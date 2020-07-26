#ifndef BARCODE_DECONVOLUTION_STAGE_HPP
#define BARCODE_DECONVOLUTION_STAGE_HPP

#include "assembly_graph/paths/path_processor.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/reads/osequencestream.hpp"
#include "pipeline/stage.hpp"
#include "sequence/sequence_tools.hpp"
#include <algorithm>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

namespace debruijn_graph {

    /*
     Minimal required definitions for the read cloud deconvolution module. 
    */

    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> lib_t;
    void processReads(conj_graph_pack &graph_pack, const lib_t &lib_10x, int i);

    class BarcodeDeconvolutionStage : public spades::AssemblyStage {
    public:
        BarcodeDeconvolutionStage() : AssemblyStage("Ariadne", "barcode_deconvolution") {}
        void run(conj_graph_pack &gp, const char*);
    };

    void BarcodeDeconvolutionStage::run(conj_graph_pack &gp, const char*)
    // Run the read cloud deconvolution module processReads() for every read library in the dataset.
    {
        gp.EnsureIndex();
        if(!gp.kmer_mapper.IsAttached()) gp.kmer_mapper.Attach();
        INFO("Read cloud deconvolution starting");
        int i = 0;
        config::dataset& dataset_info = cfg::get_writable().ds;
        for (auto& lib : dataset_info.reads) {
            if(i == 0) processReads(gp, lib, i);
            ++i;
        }
        INFO("Read cloud deconvolution finished");
    }

    /*
     Pre-processing tools to extract barcode from different read data-structures.
    */

    std::string GetTenXBarcodeFromRead(const io::PairedRead &read) {
        std::string delimiter = "BX:Z:";
        std::string delimiter2 = "-1";
        size_t start_pos = read.first().name().find(delimiter);
        size_t delimiter_size = delimiter.length();
        if (start_pos != string::npos) {
            std::string barcode = read.first().name().substr(start_pos + delimiter_size);
            size_t end_pos = barcode.find(delimiter2);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            return barcode;
        }
        return "";
    }

    std::string GetTenXBarcodeFromRead(const std::string &read) {
        std::string delimiter = "BX:Z:";
        std::string delimiter2 = "-1";
        size_t start_pos = read.find(delimiter);
        size_t delimiter_size = delimiter.length();
        if (start_pos != string::npos) {
            std::string barcode = read.substr(start_pos + delimiter_size);
            size_t end_pos = barcode.find(delimiter2);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            return barcode;
        }
        return "";
    }

    void MakeRead(io::SingleRead& read,
                  debruijn_graph::conj_graph_pack &gp,
                  std::vector<std::vector<int>>& connected_components,
                  std::vector<std::pair<int, MappingPath<EdgeId>>>& mapping_record,
                  int read_index, int complement_index)
    // Store the reads in the analysis structure (read_cloud and connected_components) and the reporting structures (*_record).
    {
        auto mapper = MapperInstance(gp);
        auto path = mapper->MapRead(read);
        mapping_record.emplace_back( read_index, path );
        connected_components.emplace_back( std::vector<int>() );
        connected_components.back().push_back( complement_index ); // Uniquely each read to the other in the pair.
    }

    /*
     Deconvolution tools applied to each cloud to subdivide them into enhanced clouds that should originate from a single 10x fragment.
    */

    std::vector<VertexId> VerticesReachedFrom(VertexId& start_vertex,
                                              debruijn_graph::conj_graph_pack &gp, int edge_size)
    // Find all vertices reachable within the search distance starting from the 3'-most vertex of the read's mapping path.
    {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(start_vertex);

        return bounded_dijkstra.ReachedVertices();
    }
    std::vector<VertexId> ConjugateVerticesReachedFrom(VertexId& start_vertex,
                                                       debruijn_graph::conj_graph_pack &gp, int edge_size)
    // Find all vertices reachable within the search distance with the 3'-most vertex of the read's mapping path as the destination.
    {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(start_vertex);

        return bounded_dijkstra.ReachedVertices();
    }

    void ClusterReads(std::vector<std::vector<int>>& connected_components,
                      std::vector<std::pair<int, MappingPath<EdgeId>>>& mapping_record,
                      debruijn_graph::conj_graph_pack &gp)
    // Generates all connections between reads in the same read cloud by constructing a Djikstra graph from accessible vertices and traversing it
    // to find other reads with mapping paths containing those vertices.
    {
        for (auto& path_1 : mapping_record) { // For each read i...
            auto read_1 = path_1.second;
            if (read_1.size()){  // Check to see if the first read has a mapping path.
                EdgeId read_1_edge_end = read_1.edge_at(read_1.size() - 1);
                EdgeId read_1_edge_start = read_1.edge_at(0);

                // Look for the set of all vertices reachable from read 1's 3'-most vertex.
                VertexId startVertex = gp.g.EdgeEnd(read_1.back().first); // 3'-most vertex of the read's mapping path.
                std::vector<VertexId> reached_vertices;
                int endDist = gp.g.length(read_1_edge_end) - read_1.end_pos(); // Distance between end of read and 3'-most vertex.
                int reducedEndDist = cfg::get().barcode_distance - endDist; // The amount of distance to search starting from the 3'-most vertex.
                if (reducedEndDist > 0){
                    reached_vertices = VerticesReachedFrom(startVertex, gp, reducedEndDist); // Find the list of vertices that can be reached with a read of MappingPath<EdgeId>*.;
                }
                std::sort(reached_vertices.begin(), reached_vertices.end()); // A sorted list of VertexIDs, which are just integers.

                // Look for the set of all vertices that can reach read 1's 5'-most vertex.
                VertexId conjStartVertex = gp.g.EdgeStart(read_1.front().first); // 5'-most vertex of the read's mapping path.
                std::vector<VertexId> conjugate_reached_vertices;
                int startDist = read_1.mapping_at(0).mapped_range.start_pos;
                int reducedStartDist = cfg::get().barcode_distance - startDist;
                if (reducedStartDist > 0){
                    conjugate_reached_vertices = ConjugateVerticesReachedFrom(conjStartVertex, gp, reducedStartDist);
                }
                std::sort(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end());

                for(auto& path_2 : mapping_record) { // ...check to see if all other reads j can be connected to it.
                    auto &read_2 = path_2.second;
                    if ( path_1.first < path_2.first && read_2.size() ){ // Omit previous reads. Check to see if the second read has a mapping path.
                        EdgeId read_2_edge_end = read_2.edge_at(read_2.size() - 1);
                        EdgeId read_2_edge_start = read_2.edge_at(0);

                        if (read_1_edge_end == read_2_edge_start ||
                            read_1_edge_end == gp.g.conjugate(read_2_edge_end)) {
                            // Does read j start on the same edge that read i ends?
                            long int distance_between_reads = read_2.start_pos() - read_1.end_pos();
                            if (std::abs(distance_between_reads) < cfg::get().barcode_distance && distance_between_reads >= 0){
                                connected_components[ path_1.first ].push_back( path_2.first );
                                connected_components[ path_2.first ].push_back( path_1.first );
                            }
                        } else if (read_1_edge_start == read_2_edge_end ||
                                   read_1_edge_start == gp.g.conjugate(read_2_edge_start) ) {
                            // Does read j end on the same edge that read i starts?
                            int path2_conj_end = gp.g.length(read_2.edge_at(read_2.size() - 1)) - read_2.end_pos();
                            long int distance_between_reads = path2_conj_end - read_1.end_pos();
                            if (std::abs(distance_between_reads) < cfg::get().barcode_distance && distance_between_reads >= 0){
                                connected_components[ path_1.first ].push_back( path_2.first );
                                connected_components[ path_2.first ].push_back( path_1.first );
                            }
                        }
                        else { // Otherwise, is the read between read j traversable within the Djikstra graph?
                            // Evaluate the 5' and 3'-most vertices of read 2.
                            VertexId forward_end_three = gp.g.EdgeEnd(read_2_edge_end);
                            // Assuming read is same-strand. The 3'-most vertex of the last edge that read 2 sits on.
                            VertexId reverse_end_three = gp.g.conjugate(gp.g.EdgeStart(read_2_edge_start));
                            // Assuming read is opposite-strand. The 5'-most vertex of the first edge that read 2 sits on.
                            VertexId forward_end_five = gp.g.EdgeStart(read_2_edge_start);
                            // Assuming read is same-strand. The 5'-most vertex of the last edge that read 2 sits on.
                            VertexId reverse_end_five = gp.g.conjugate(gp.g.EdgeEnd(read_2_edge_end));
                            // Assuming read is opposite-strand. The 3'-most vertex of the last edge that read 2 sits on.
                            if (std::binary_search(reached_vertices.begin(), reached_vertices.end(), forward_end_five) ||
                                std::binary_search(reached_vertices.begin(), reached_vertices.end(), reverse_end_five)) {
                                // If the 5'-most vertex can be found in the Djikstra graph, then add an edge between the two reads.
                                connected_components[ path_1.first ].push_back( path_2.first );
                                connected_components[ path_2.first ].push_back( path_1.first );
                            } else if (std::binary_search(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end(), forward_end_three) ||
                                       std::binary_search(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end(), forward_end_three)) {
                                // If the 3'-most vertex can be found in the conjugate Djikstra graph, then add an edge between the two reads.
                                connected_components[ path_1.first ].push_back( path_2.first );
                                connected_components[ path_2.first ].push_back( path_1.first );
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     Reporting tools to output enhanced read clouds.
    */

    inline std::string GetUpdatedReadName(std::string& read, size_t cloud_number)
    // Updates original read cloud number to enhanced read cloud number.
    {
        std::string delimiter = "-1";
        size_t start_pos = read.find(delimiter);
        size_t delimiter_size = delimiter.length();
        if (start_pos != string::npos) {
            std::string read_name = read.substr(0, start_pos) + "-" + std::to_string(cloud_number);
//            if (read.substr(start_pos+delimiter_size).find("_RC") != std::string::npos) { // If it is the second read, tag with "_RC".
//                read_name.append("_RC");
//            } // Got rid of these because the paired reads are de-interleaved in the output step.
            return read_name;
        }
        return "";
    }

    void OutputPaths(std::vector<std::vector<int>>& connected_components,
                     std::vector<std::pair<int, MappingPath<EdgeId>>>& mapping_record,
                     std::vector<std::tuple<std::string, std::string, std::string>>& read_record,
                     std::string& barcode,
                     io::OFastqReadStream& fastq_stream_forward,
                     io::OFastqReadStream& fastq_stream_reverse,
                     std::ofstream& enh_cloud_stats)
    // Report enhanced read clouds, which are depth-first traversals through sets of connected reads.
    {
        std::vector<int> cloud_queue;
        std::vector<int> visited_reads;
        std::vector<std::tuple<std::string, std::string, std::string>> enhanced_record;
        std::vector<std::tuple<std::string, std::string, std::string>> master_record;
        size_t num_enhanced_cloud = 1;

        for (size_t i = 0; i < mapping_record.size(); ++i) {
            auto read_find = std::find(visited_reads.begin(), visited_reads.end(), mapping_record[i].first);
            if ( read_find != visited_reads.end() ) { // If this read has already been grouped, skip.
                continue;
            }
            else {
                size_t num_reads = 0; // The number of reads in the enhanced cloud.
                visited_reads.push_back( mapping_record[i].first ); // Add read to vector of already-grouped reads.
                cloud_queue.push_back( mapping_record[i].first );
                std::vector<std::tuple<std::string, std::string, std::string>> temporary_enhanced;
                do {
                    int& cloud_index = cloud_queue.front(); // Get the unique index of read j.
                    temporary_enhanced.emplace_back( read_record[cloud_index] );
                    ++num_reads;
                    for ( int j = 0; j < connected_components[cloud_index].size(); ++j ) // Do this for every read k that read j is connected to.
                    {
                        auto current_index = std::find(visited_reads.begin(), visited_reads.end(), connected_components[cloud_index][j]);
                        if( current_index == visited_reads.end() ){ // If the current read has not been grouped...
                            auto additional_read = connected_components[cloud_index][j];
                            cloud_queue.emplace_back( additional_read ); // ...add it to the cloud queue
                            visited_reads.emplace_back( additional_read ); // ...add it to the vector of grouped reads.
                        }
                    }
                    cloud_queue.erase(cloud_queue.begin());
                }
                while ( !cloud_queue.empty() );
                if ( num_reads > 2 ) { // This enhanced read cloud contains more than just a pair of reads.
                    for ( auto& read : temporary_enhanced ){
                        enhanced_record.emplace_back( std::make_tuple(
                                GetUpdatedReadName( std::get<0>(read), num_enhanced_cloud ),
                                std::get<1>(read), std::get<2>(read) ) );
                    }
                    enh_cloud_stats << barcode << "," << num_enhanced_cloud << "," << num_reads << endl; // Reporting enhanced read cloud information.
                    ++num_enhanced_cloud;
                } else { // Just a pair of reads. Add to the un-enhanced cluster of reads.
                    for ( auto& read : temporary_enhanced ){
                        master_record.emplace_back( std::make_tuple(
                                GetUpdatedReadName( std::get<0>(read), 0),
                                std::get<1>(read), std::get<2>(read) ) );
                    }
                }
            }
        }
        enh_cloud_stats << barcode << "," << 0 << "," << master_record.size() << endl; // Reporting un-enhanced read cloud information.

        // Add reads in the enhanced read clouds to the master collection of reads.
        while ( !enhanced_record.empty() ){
            master_record.emplace_back( enhanced_record.front() );
            enhanced_record.erase(enhanced_record.begin());
        }

        int direction = 0;
        while ( !master_record.empty() ){
            auto read = master_record.front();
            io::SingleRead input(std::get<0>(read), std::get<1>(read), std::get<2>(read));
            if ( direction == 0 || direction % 2 == 0 ) { // If first read in pair, put in forward file.
                fastq_stream_forward << input;
            } else { // If second read in pair, put in reverse complement file.
                fastq_stream_reverse << input;
            }
            ++direction;
            connected_components.erase(connected_components.begin());
            mapping_record.erase(mapping_record.begin());
            read_record.erase(read_record.begin());
            master_record.erase(master_record.begin());
        }
    }

    /*
     main() function.
    */

    void processReads(debruijn_graph::conj_graph_pack &graph_pack, const lib_t& lib_10x, const int i)
    // Main coordinating function. Separates reads into original read clouds and maps them to the assembly graph for deconvolution.
    {
        // Tools to locate reads along the compacted de Bruijn graph.
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;

        // Setting the output files.
        std::string file_name = cfg::get().output_dir + std::to_string(cfg::get().barcode_distance)  + "_" + std::to_string(i) ;
        io::OFastqReadStream fastq_stream_forward(file_name + ".R1.fastq");
        io::OFastqReadStream fastq_stream_reverse(file_name + ".R2.fastq");
        std::ofstream enh_cloud_stats;
        enh_cloud_stats.open(file_name + ".summary.txt", std::ofstream::out | std::ofstream::app);

        // Collection of read information to be outputted, in the original read cloud, and enhanced read clouds respectively.
        std::vector<std::vector<int>> connected_components;
        std::vector<std::pair<int, MappingPath<EdgeId>>> mapping_record;
        std::vector<std::tuple<std::string, std::string, std::string>> read_record;
        std::vector<std::tuple<std::string, std::string, std::string>> unbarcoded_record;

        int num_original_clouds = 1;
        int num_cloud_reads = 0;
        std::string current_barcode = "";

        while ( !stream->eof() ) {
            *stream >> read;
            std::string barcode_string = GetTenXBarcodeFromRead(read);
            if ( !barcode_string.empty() ){
                if( barcode_string != current_barcode && !mapping_record.empty() ){ // This read belongs to the next read-cloud, process the current read cloud.
                    num_original_clouds++;
                    if( current_barcode.empty() ){ // For the first original read cloud.
                        current_barcode = GetTenXBarcodeFromRead(std::get<0>(read_record[0]));
                    }
                    if(num_original_clouds % 10000 == 0 ) {
                        INFO(num_original_clouds << ": Processing read cloud " << current_barcode << ", " << mapping_record.size() << " reads");
                    }
                    ClusterReads(connected_components, mapping_record, graph_pack);
                    OutputPaths(connected_components, mapping_record, read_record, current_barcode, fastq_stream_forward, fastq_stream_reverse, enh_cloud_stats);
                    num_cloud_reads = 0;
                }
                // For each pair of reads, map them to the assembly graph and extract their FastQ information.
                read_record.emplace_back( std::make_tuple( read.first().name(), read.first().GetSequenceString(),
                                                           read.first().GetQualityString() ) );
                read_record.emplace_back( std::make_tuple( read.second().name(), Complement(read.second().GetSequenceString()),
                                                           Reverse(read.second().GetQualityString()) ) );
                MakeRead(read.first(), graph_pack, connected_components, mapping_record, num_cloud_reads, num_cloud_reads + 1);
                MakeRead(read.second(), graph_pack, connected_components, mapping_record, num_cloud_reads + 1, num_cloud_reads);

                num_cloud_reads = num_cloud_reads + 2;
                current_barcode = barcode_string;
            }
            else { // This read does not have a barcode, and will be reported at the end.
                unbarcoded_record.emplace_back( std::make_tuple( read.first().name(), read.first().GetSequenceString(),
                                                                 read.first().GetQualityString() ) );
                unbarcoded_record.emplace_back( std::make_tuple( read.second().name(), Complement(read.second().GetSequenceString()),
                                                                 Reverse(read.second().GetQualityString()) ) );
            }
        }
        INFO(num_original_clouds << " original read clouds in total");
        ClusterReads(connected_components, mapping_record, graph_pack);
        OutputPaths(connected_components, mapping_record, read_record, current_barcode, fastq_stream_forward, fastq_stream_reverse, enh_cloud_stats);

        int direction = 0;
        while ( !unbarcoded_record.empty() ){ // These reads do not have a barcode, and are appended to the end of the FastQ.
            auto read = unbarcoded_record.front();
            io::SingleRead input(std::get<0>(read), std::get<1>(read), std::get<2>(read));
            if ( direction == 0 || direction % 2 == 0 ) { // If first read in pair, put in forward file.
                fastq_stream_forward << input;
            } else { // If second read in pair, put in reverse complement file.
                fastq_stream_reverse << input;
            }
            ++direction;
            unbarcoded_record.erase(unbarcoded_record.begin());
        }
        enh_cloud_stats << "no_barcode,NA," << direction << endl; // Reporting enhanced read cloud information.
        enh_cloud_stats.close();
    }
}

#endif