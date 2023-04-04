#ifndef BARCODE_DECONVOLUTION_STAGE_HPP
#define BARCODE_DECONVOLUTION_STAGE_HPP

#include "assembly_graph/paths/path_processor.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/osequencestream.hpp"
#include "pipeline/stage.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/memory_limit.hpp"
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <numeric>
#include <omp.h>
#include <string>
#include <ostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <vector>

namespace debruijn_graph {

    typedef std::tuple<std::string, std::string, std::string> fq_info;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> lib_t;

    class BarcodeDeconvolutionStage : public spades::AssemblyStage {
    public:
        BarcodeDeconvolutionStage() : AssemblyStage("Ariadne", "barcode_deconvolution") {}
        void run(conj_graph_pack &gp, const char*);
    };

    // Step -1A: Set up the output streams
    std::ofstream MakeOutputStream(std::string suffix) {
        std::ofstream stream;
        const size_t bufsize = 1024*1024;
        char buf[bufsize];
        stream.rdbuf()->pubsetbuf(buf, bufsize);
        std::string prefix = cfg::get().output_dir + std::to_string(cfg::get().search_distance);
        stream.open(prefix + suffix, std::ofstream::out);
        return stream;
    }

    // Step -1B: Count number of reads in the dataset.
    int CountReads(const lib_t& lib_10x)
    {
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;
        int num_pairs_total = 0;

        while ( !stream->eof() ) {
            *stream >> read;
            num_pairs_total++;
        }
        stream->close();
        INFO(num_pairs_total << " read-pairs to process");
        return num_pairs_total;
    }

    // Step -1C: Get barcode (currently, only 10X format supported)
    inline std::string GetBarcode(std::string &read) 
    {
        std::string start = "BX:Z:";
        std::string end = "-1";
        size_t start_pos = read.find(start);
        size_t delimiter_size = start.length();
        if (start_pos != string::npos) {
            std::string barcode = read.substr(start_pos + delimiter_size);
            size_t end_pos = barcode.find(end);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            return barcode;
        }
        return "";
    }

    // Step 0A: Get barcode (currently, only 10X format supported)
    inline std::string GetBarcode(const io::PairedRead &read) 
    {
        std::string start = "BX:Z:";
        std::string end = "-1";
        size_t start_pos = read.first().name().find(start);
        size_t delimiter_size = start.length();
        if (start_pos != string::npos) {
            std::string barcode = read.first().name().substr(start_pos + delimiter_size);
            size_t end_pos = barcode.find(end);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            return barcode;
        }
        return "";
    }

    // Step 0B: Store the reads in the analysis structures (reads and edges) and the reporting structure (infos).
    void MakeRead(io::SingleRead& read,
                  debruijn_graph::conj_graph_pack& gp,
                  std::vector<std::vector<MappingPath<EdgeId>>>& subsection_edges,
                  std::vector<std::vector<fq_info>>& subsection_infos,
                  fq_info& read_infos, 
                  int read_idx)
    {
        auto mapper = MapperInstance(gp);
        auto path = mapper->MapRead(read);
        // INFO(std::get<0>(read_infos) << " " << path);
        subsection_edges.back().emplace_back( path );
        subsection_infos.back().emplace_back( read_infos ); 
    }

    // Step 0: Load reads in PairedEasyStream into read cloud data structures
    int LoadReads(std::vector<std::vector<MappingPath<EdgeId>>>& subsection_edges,
                  std::vector<std::vector<fq_info>>& subsection_infos,
                  debruijn_graph::conj_graph_pack& gp, const lib_t& lib_10x, 
                  int pair_start_idx, int num_loadable_clouds)
    {
        // Tools to locate reads along the compacted de Bruijn graph.
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;

        // Initialize barcode string and cloud-tracker values.
        std::string current_barcode = "";
        int pair_idx = 0;
        int num_clouds = 0;
        int num_pairs = 0; 

        while ( !stream->eof() ) {
            *stream >> read;
            std::string barcode = GetBarcode(read);
            if ( pair_idx < pair_start_idx || barcode.empty() ) {                   // If this read has already been loaded or doesn't have barcodes 
                pair_idx++;
                continue;
            }
            if ( barcode != current_barcode ){                                      // This barcoded read belongs to the next read-cloud. Start new storage objects.
                if ( num_clouds == num_loadable_clouds ) {                          // If no more loadable clouds, stop.
                    break;
                }
                num_clouds++;
                num_pairs = 0; 
                subsection_edges.emplace_back( std::vector<MappingPath<EdgeId>>() );
                subsection_infos.emplace_back( std::vector<fq_info>() );
            }
            // For each pair of reads, map them to the assembly graph and extract their FastQ information.
            auto fwd_read = read.first();
            fq_info fwd_info = std::make_tuple( fwd_read.name(), fwd_read.GetSequenceString(), fwd_read.GetPhredQualityString() );
            int fwd_idx = num_pairs * 2;
            MakeRead(fwd_read, gp, subsection_edges, subsection_infos, fwd_info, fwd_idx );
            auto rev_read = read.second();
            fq_info rev_info = std::make_tuple( rev_read.name(), Reverse(Complement(rev_read.GetSequenceString())), Reverse(rev_read.GetPhredQualityString()) );
            int rev_idx = num_pairs * 2 + 1;
            MakeRead(rev_read, gp, subsection_edges, subsection_infos, rev_info, rev_idx );

            current_barcode = barcode;
            pair_idx++;
            num_pairs++;
        }
        stream->close();
        return pair_idx;
    }

    // Step 1Ai: Find all vertices reachable within the search distance going 5' from a vertex.
    std::vector<VertexId> DjikstraSearch(VertexId& vertex_start,
                                         debruijn_graph::conj_graph_pack &gp, int edge_size)
    {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(vertex_start);

        return bounded_dijkstra.ReachedVertices();
    }

    // Step 1Aii: Find all vertices reachable within the search distance going 3' from a vertex.
    std::vector<VertexId> DjikstraConjSearch(VertexId& vertex_start,
                                             debruijn_graph::conj_graph_pack &gp, int edge_size)
    {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(vertex_start);

        return bounded_dijkstra.ReachedVertices();
    }

    // Step 1A: Find all vertices reachable within the search distance from a read's mapping position.
    std::vector<VertexId> GetNeighbors(MappingPath<EdgeId>& edge, debruijn_graph::conj_graph_pack &gp)
    {
        int search_dist = cfg::get().search_distance;
        std::vector<VertexId> vertex_vct; 

        if ( ! edge.empty() ) {
            EdgeId edge_end = edge.edge_at(edge.size() - 1);
            VertexId vertex_five = gp.g.EdgeStart(edge.front().first);              // 5'-most vertex of the read's mapping path.
            // INFO(edge_end << " " << vertex_five);

            vertex_vct.emplace_back(vertex_five);                                   // Add the read's 5'- and 3'-most vertices to the set

            VertexId vertex_three = gp.g.EdgeEnd(edge.back().first);                // 3'-most vertex of the read's mapping path
            vertex_vct.emplace_back(vertex_three);
            int dist_to_three = gp.g.length(edge_end) - edge.end_pos();             // Distance between end of read and 3'-most vertex.
            int reduced_fwd_dist = search_dist - dist_to_three;                     // The amount of distance to search starting from the 3'-most vertex.
            if (reduced_fwd_dist > 0){
                std::vector<VertexId> vertex_tmp = DjikstraSearch(vertex_three, gp, reduced_fwd_dist);        // Find the set of vertices that can be reached with a read on MappingPath<EdgeId>*
                vertex_vct.insert( vertex_vct.end(), vertex_tmp.begin(), vertex_tmp.end() );
                // INFO("Forward " << vertex_tmp);
            }
            int dist_to_five = edge.mapping_at(0).mapped_range.start_pos;           // Distance between start of read and 5'-most vertex.
            int reduced_rev_dist = search_dist - dist_to_five;             
            if (reduced_rev_dist > 0){
                std::vector<VertexId> vertex_tmp = DjikstraConjSearch(vertex_five, gp, reduced_rev_dist);     // Find the set of vertices that can be reached with a read on MappingPath<EdgeId>*
                vertex_vct.insert( vertex_vct.end(), vertex_tmp.begin(), vertex_tmp.end() );
            }
        }
        return vertex_vct;
    }

    // Step 1B: See if there is overlap between the new read's set of reachable vertices and the current subcloud's list of vertices
    bool HasIntersection(std::vector<VertexId>& read_vct, std::vector<VertexId>& subcloud_vct)
    {
        for (size_t i = 0; i < read_vct.size(); ++i) {
            if (std::binary_search(subcloud_vct.begin(), subcloud_vct.end(), read_vct[i])) {
                return true;
            }
        }
        return false;
    }

    // Step 1: Find overlaps between the sets of vertices that the reads
    std::vector<std::vector<int>> MakeClouds(std::vector<MappingPath<EdgeId>>& edges,
                                             debruijn_graph::conj_graph_pack &gp)
    {
        std::vector<std::vector<VertexId>> subcloud_vert_lst;                       // subcloud[subcloud[vertex_id]]
        std::vector<std::vector<int>> raw_subcloud_lst;                             // subcloud[subcloud{read_id]]
        std::vector<int> unmappable_read_lst;                                       // subcloud[read_id]. Keeps unmappable reads out of the to-check list
        
        for (size_t i = 0; i < edges.size(); i+=2) {                                // For each pair of reads with indices [i, i + 1]...
            // INFO("Reads " << i << " and " << i+1);
            std::vector<VertexId> vertex_fwd = GetNeighbors(edges[i], gp); 
            std::vector<VertexId> vertex_rev = GetNeighbors(edges[i + 1], gp); 
            std::vector<VertexId> vertex_vct;                                       // Initialize an empty vertex set

            vertex_vct.insert( vertex_vct.end(), vertex_fwd.begin(), vertex_fwd.end() );
            vertex_vct.insert( vertex_vct.end(), vertex_rev.begin(), vertex_rev.end() );
            std::sort(vertex_vct.begin(), vertex_vct.end());
            // INFO(vertex_vct);

            std::vector<int> read_id_vct = {i, i + 1};
            if ( ! vertex_vct.empty() ) {                                              // Check to see if the pair has a mapping path.
                std::vector<int> cloud_intersections; 
                for (size_t j = 0; j < subcloud_vert_lst.size(); ++j) {             // For each subcloud vertex set index j...
                    if ( ! subcloud_vert_lst[j].empty() && HasIntersection(vertex_vct, subcloud_vert_lst[j]) ){       // If the read vertices and subcloud vertices has overlap
                        cloud_intersections.emplace_back( j );
                    }
                }
                if ( cloud_intersections.size() == 0 ) {                                // The read vertices and all subcloud vertices have no overlap
                    subcloud_vert_lst.emplace_back( vertex_vct );                       // Initialize the sets of vertices and read ID as a new subcloud
                    raw_subcloud_lst.emplace_back( read_id_vct );
                    // INFO("New cloud");
                } else  {                                                               // The read vertices overlap with at least one set of subclouds
                    int first_idx = cloud_intersections[0];
                    // INFO(edges.size() << " " << i << " " << subcloud_vert_lst[first_idx].size() << " " << subcloud_vert_lst[first_idx]);
                    // INFO(edges.size() << " " << i << " " << raw_subcloud_lst[first_idx].size() << " " << raw_subcloud_lst[first_idx]);
                    if ( cloud_intersections.size() > 1 ) {                                    // The read vertices overlap with multiple sets of subcloud vertices- join them!
                        for (size_t j = 1; j < cloud_intersections.size(); ++j) { 
                            int curr_idx = cloud_intersections[j];
                            // INFO(edges.size() << " " << i << " " << curr_idx);
                            subcloud_vert_lst[first_idx].insert( subcloud_vert_lst[first_idx].end(), subcloud_vert_lst[curr_idx].begin(), subcloud_vert_lst[curr_idx].end() );
                            raw_subcloud_lst[first_idx].insert( raw_subcloud_lst[first_idx].end(), raw_subcloud_lst[curr_idx].begin(), raw_subcloud_lst[curr_idx].end() );
                            std::vector<VertexId> tmp_vct;
                            subcloud_vert_lst[curr_idx] = tmp_vct;
                            std::vector<int> tmp_reads; 
                            raw_subcloud_lst[curr_idx] = tmp_reads;
                        }
                    }
                    subcloud_vert_lst[first_idx].insert( subcloud_vert_lst[first_idx].end(), vertex_vct.begin(), vertex_vct.end() );
                    std::sort(subcloud_vert_lst[first_idx].begin(), subcloud_vert_lst[first_idx].end());
                    raw_subcloud_lst[first_idx].insert( raw_subcloud_lst[first_idx].end(), read_id_vct.begin(), read_id_vct.end() ); 
                    // INFO(edges.size() << " " << i << " " << subcloud_vert_lst[first_idx].size() << " " << subcloud_vert_lst[first_idx]);
                    // INFO(edges.size() << " " << i << " " << raw_subcloud_lst[first_idx].size() << " " << raw_subcloud_lst[first_idx]);
                }
            } else {
                unmappable_read_lst.insert( unmappable_read_lst.end(), read_id_vct.begin(), read_id_vct.end() ); 
            }
            // INFO(i << " " << i + 1);
            // INFO(subcloud_vert_lst);
            // INFO(raw_subcloud_lst);
            // INFO(unmappable_read_lst);
        }

        // Add subclouds that are larger in pairs to the final report list, and the orphan pairs to the unmappable list
        std::vector<std::vector<int>> subcloud_read_lst; 
        for (size_t i = 0; i < raw_subcloud_lst.size(); ++i) {                      // For each subcloud index i in reverse order...
            if ( raw_subcloud_lst[i].empty() ){
                // INFO(edges.size() << " " << i << " empty set");
                continue;
            }
            else if ( raw_subcloud_lst[i].size() == 2 ) {                               // If this is a disconnected read-pair...
                // INFO(edges.size() << " " << i << " orphaned");
                unmappable_read_lst.insert( unmappable_read_lst.end(), raw_subcloud_lst[i].begin(), raw_subcloud_lst[i].end() ); 
                // INFO(i << " orphan " << raw_subcloud_lst[i]);
            } else {
                // INFO(edges.size() << " " << i << " " << subcloud_read_lst.size() << " " << raw_subcloud_lst[i].size());
                subcloud_read_lst.emplace_back( raw_subcloud_lst[i] );
                // INFO(i << " new subcloud " << raw_subcloud_lst[i]);
            }
        }
        subcloud_read_lst.insert(subcloud_read_lst.begin(), unmappable_read_lst);   // Unmappable/orphan reads at the beginning 
        // INFO(unmappable_read_lst);
        // INFO(subcloud_read_lst);

        // INFO(raw_subcloud_lst);
        return subcloud_read_lst;
    }

    // Step 3A: Concatenates read information into FastQ format
    inline std::string MakeFastQString(std::string& name, std::string& seq, std::string& qual)
    {
        std::string read = "@" + name + "\n" + seq + "\n+\n" + qual + "\n";
        return read;
    }

    // Step 3: Update each read's barcode with the subcloud ID and make the 
    std::vector<std::string> UpdateReads(std::vector<fq_info>& infos,
                     std::vector<std::vector<int>>& raw_subcloud_lst)              // Note: This already has reads sorted in barcode order!
    {
        std::string fwd_str = "";
        std::string rev_str = "";
        size_t start_idx = 0; 
        if ( raw_subcloud_lst[0].empty() ) {                                       // If there are no unmapped and/or disconnected reads...
            start_idx = 1;
        }
        for (size_t i = start_idx; i < raw_subcloud_lst.size(); ++i) {             // For each subcloud index i...
            std::string subcloud_idx = std::to_string(i);
            for (size_t j = 0; j < raw_subcloud_lst[i].size(); ++j) {              // For each read index at j...
                auto old_info = infos[raw_subcloud_lst[i][j]];
                std::string old_name = std::get<0>(old_info);
                std::string new_name = old_name; 
                size_t start_pos = old_name.find("-1");                             // Position in the read where the subcloud ID starts
                if ( start_pos != std::string::npos ) {                             // If the read has a barcode, replace '-1' with the subcloud ID
                    new_name = old_name.substr(0, start_pos) + "-" + subcloud_idx;
                }
                std::string new_info = MakeFastQString(new_name, std::get<1>(old_info), std::get<2>(old_info)); 
                if ( j % 2 == 0 ) {                                                     // If the read is the first in a pair...
                    fwd_str.append(new_info);
                    // INFO(new_info);
                } else {
                    rev_str.append(new_info);
                }
            }
        }
        std::vector<std::string> fq_output{ fwd_str, rev_str };
        return fq_output;
    }


    // Step 4: Write enhanced reads to FastQs
    void OutputReads(std::vector<std::string>& fq_output,
                     std::ofstream& fwd_fq_stream,
                     std::ofstream& rev_fq_stream)
    {
        fwd_fq_stream << fq_output[0] << std::flush;
        rev_fq_stream << fq_output[1] << std::flush;
    }

    /* run() */
    void BarcodeDeconvolutionStage::run(conj_graph_pack &gp, const char*)
    {
        gp.EnsureIndex();
        if (!gp.kmer_mapper.IsAttached()) gp.kmer_mapper.Attach();
        INFO("Read cloud deconvolution starting");
        config::dataset& dataset_info = cfg::get_writable().ds;
        lib_t& lib_10x = dataset_info.reads[0];

        // Output stream objects.
        std::ios::sync_with_stdio(false);
        std::cin.tie(nullptr);
        std::ofstream fwd_fq_stream = MakeOutputStream(".R1.fastq");
        std::ofstream rev_fq_stream = MakeOutputStream(".R2.fastq");

        // Collection of read information from the original read clouds
        std::vector<std::vector<MappingPath<EdgeId>>> subsection_edges;
        std::vector<std::vector<fq_info>> subsection_infos;

        // Thread- and read-based limits on the loop
        int pair_start_idx = 0;
        int num_pairs_total = CountReads(lib_10x);                                  // Count the number of reads in the dataset
        int num_loadable_clouds = (int)( (utils::get_free_memory() * 1.0) / 150000.0 ); // Available memory in bytes / bytes/cloud
        int cloud_size_filter = cfg::get().size_cutoff;
        INFO("Cutoff for original cloud size is " << cloud_size_filter << " reads");

        while ( pair_start_idx < num_pairs_total ) {
            INFO("Loading original read clouds from library");
            int pair_end_idx = LoadReads(subsection_edges, subsection_infos, gp, lib_10x, pair_start_idx, num_loadable_clouds);
            size_t num_clouds = subsection_edges.size();
            INFO(num_clouds << " read clouds loaded");
#pragma omp parallel for shared(subsection_edges, subsection_infos) schedule(dynamic, 1) num_threads(cfg::get().max_threads)
            for (size_t c = 0; c < num_clouds; ++c) {
                std::string barcode = GetBarcode(std::get<0>(subsection_infos[c][0]));
                int num_reads = subsection_edges[c].size(); 
                std::vector<std::vector<int>> subcloud_read_lst;
                if ( ! barcode.empty() && num_reads >= cloud_size_filter ) {
                    if (c % 50000 == 0) {
                        INFO(c << ": Ariadne deconvolving " << subsection_edges[c].size() << " reads on thread " << omp_get_thread_num());
                    }
                    subcloud_read_lst = MakeClouds(subsection_edges[c], gp);
                } else {
                    std::vector<int> tmp(num_reads, 0);                             // per_read[subcloud_id], default is 0
                    subcloud_read_lst.emplace_back( tmp );
                }
                std::vector<std::string> fq_output = UpdateReads(subsection_infos[c], subcloud_read_lst);
#pragma omp critical
                {
                    OutputReads(fq_output, fwd_fq_stream, rev_fq_stream);
                }
                // if (c % 1 == 0) {
                //     INFO(c << ": Ariadne finished " << subsection_edges[c].size() << " reads on thread " << omp_get_thread_num());
                // }
            }
            subsection_edges.clear();
            subsection_infos.clear();
            // INFO(pair_start_idx << " " << pair_end_idx);
            pair_start_idx = pair_end_idx;                                         // Start after this read index next time!
        }
        fwd_fq_stream.close();
        rev_fq_stream.close();
        INFO("Read cloud deconvolution finished");
    }

}

#endif
