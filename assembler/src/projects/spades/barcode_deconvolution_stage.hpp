#ifndef BARCODE_DECONVOLUTION_STAGE_HPP
#define BARCODE_DECONVOLUTION_STAGE_HPP

#include "pipeline/stage.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include <unordered_map>
#include <string>
#include <cstdlib>
// LM: All changes introduced by Lauren Mak will be commented and tagged with LM.

namespace debruijn_graph {

    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> lib_t;
    typedef std::vector<io::SequencingLibrary<debruijn_graph::config::LibraryData>> lib_vector_t;

    class SimplePathExtractor;
    class LongReadsCreator;
    std::string GetTenXBarcodeFromRead(const io::PairedRead &read);
    void processReads(debruijn_graph::conj_graph_pack &graph_pack, const lib_t& lib_10x);

    class BarcodeDeconvolutionStage : public spades::AssemblyStage {
    public:

        BarcodeDeconvolutionStage(bool preliminary = false) 
        : AssemblyStage(preliminary ? "Preliminary Barcode Deconvolution" : "BarcodeDeconvolution",
                        preliminary ? "barcode_deconvolution_preliminary" : "barcode_deconvolution") {}

        void run(conj_graph_pack &gp, const char*);

    };

    void BarcodeDeconvolutionStage::run(conj_graph_pack &gp, const char*){
        gp.EnsureIndex();
        if(!gp.kmer_mapper.IsAttached()) gp.kmer_mapper.Attach();
        INFO("Barcode Deconvolution Started");
        int i = 0;
        config::dataset& dataset_info = cfg::get_writable().ds;
        for (auto& lib : dataset_info.reads) {

            if(i == 0) processReads(gp, lib);
            ++i;
        }
        INFO("Barcode Deconvolution Ended");
    }

    std::string GetTenXBarcodeFromRead(const io::PairedRead &read) {
        std::string delimeter = "BX:Z:";
        std::string delimeter2 = "-1";
        size_t start_pos = read.first().name().find(delimeter);
        size_t delimeter_size = delimeter.length();
        if (start_pos != string::npos) {
            std::string barcode = read.first().name().substr(start_pos + delimeter_size);
            size_t end_pos = barcode.find(delimeter2);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            // TRACE(barcode);
            return barcode;
        }
        return "";
    }

    std::string GetTenXBarcodeFromRead(const std::string &read) {
        std::string delimeter = "BX:Z:";
        std::string delimeter2 = "-1";
        size_t start_pos = read.find(delimeter);
        size_t delimeter_size = delimeter.length();
        if (start_pos != string::npos) {
            std::string barcode = read.substr(start_pos + delimeter_size);
            size_t end_pos = barcode.find(delimeter2);
            if (end_pos != string::npos) {
                barcode = barcode.substr(0, end_pos);
            }
            // TRACE(barcode);
            return barcode;
        }
        return "";
    }

    void AddEdge(std::unordered_map<const MappingPath<EdgeId>*, path_extend::BidirectionalPath*>& visited, 
        const std::pair<MappingPath<EdgeId>, std::pair<std::pair<std::string, std::string>, std::string>>& path, 
        const std::pair<MappingPath<EdgeId>, std::pair<std::pair<std::string, std::string>, std::string>>& path2,
        std::unordered_map<std::string, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>>& path_set,
        debruijn_graph::conj_graph_pack &gp, 
        std::string &barcode,
        bool first,
        std::queue<path_extend::BidirectionalPath*>& stability_queue){
        if(!visited.count(&path.first)){
            std::vector<EdgeId> edges = path.first.simple_path();
            path_extend::BidirectionalPath* bidirectional_path = new path_extend::BidirectionalPath(gp.g);
            bidirectional_path->SetConjPath(new path_extend::BidirectionalPath(gp.g));
            for (auto e : edges) {
                    bidirectional_path->PushBack(e);
                    bidirectional_path->GetConjPath()->PushBack(gp.g.conjugate(e));
            }
            // TODO -stop using pairs and store header(name), quality and sequence smarter
            bidirectional_path->name = path.second.first.first;
            bidirectional_path->quality_string = path.second.first.second;
            bidirectional_path->sequence_string = path.second.second;
            //                                          Pointer to path        vector of paths with edges to first(adjacencies)
            path_set[barcode][bidirectional_path] =  std::vector<path_extend::BidirectionalPath*>();
            stability_queue.push(bidirectional_path);
            visited[&path.first] = bidirectional_path;
            first = false;
        }
        // // If Read hasn't been traversed before & corresponding bidirectional path has not been made
        // if(!visited.count(&path.first) && first){
        //     std::vector<EdgeId> edges = path.first.simple_path();
        //     path_extend::BidirectionalPath* bidirectional_path = new path_extend::BidirectionalPath(gp.g);
        //     bidirectional_path->SetConjPath(new path_extend::BidirectionalPath(gp.g));
        //     for (auto e : edges) {
        //             bidirectional_path->PushBack(e);
        //             bidirectional_path->GetConjPath()->PushBack(gp.g.conjugate(e));
        //     }
        //     // TODO -stop using pairs and store header(name), quality and sequence smarter
        //     bidirectional_path->name = path.second.first.first;
        //     bidirectional_path->quality_string = path.second.first.second;
        //     bidirectional_path->sequence_string = path.second.second;
        //     //                                          Pointer to path        vector of paths with edges to first(adjacencies)
        //     path_set[barcode][bidirectional_path] =  std::vector<path_extend::BidirectionalPath*>();
        //     stability_queue.push(bidirectional_path);
        //     visited[&path.first] = bidirectional_path;
        //     first = false;
        // }
        // } else if (first){
        //     INFO("2");
        //     //                                          Pointer to path        vector of paths with edges to first
        //     path_set[barcode][visited[&path.first]] = std::vector<path_extend::BidirectionalPath*>();
        //     first = false;
        // }
        if(!visited.count(&path2.first)){
            std::vector<EdgeId> edges = path2.first.simple_path();
            path_extend::BidirectionalPath* bidirectional_path2 = new path_extend::BidirectionalPath(gp.g);
            bidirectional_path2->SetConjPath(new path_extend::BidirectionalPath(gp.g));
            for (auto e : edges) {
                    bidirectional_path2->PushBack(e);
                    bidirectional_path2->GetConjPath()->PushBack(gp.g.conjugate(e));
            }
            bidirectional_path2->name = path2.second.first.first;
            bidirectional_path2->quality_string = path2.second.first.second;
            bidirectional_path2->sequence_string = path2.second.second;
            path_set[barcode][visited[&path.first]].push_back(bidirectional_path2);
            visited[&path2.first] = bidirectional_path2;
            stability_queue.push(bidirectional_path2);
        } else{
            // add it to the adjacency list rather than creating a new adjacency list
            path_set[barcode][visited[&path.first]].push_back(visited[&path2.first]);
        }
    }

    // void clusterReads(debruijn_graph::conj_graph_pack &gp,
    //     std::vector<std::pair<MappingPath<EdgeId>, std::pair<std::pair<std::string, std::string>, std::string>>>& paths,
    //     std::unordered_map<std::string, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>>& path_set,
    //     std::string& barcode){
    //     std::unordered_map<const MappingPath<EdgeId>*, path_extend::BidirectionalPath*> visited;
    //     for (auto const& path : paths) {
    //         bool first = true;
    //         for(auto const& path2 : paths) {

    //             if(&path.first != &path2.first){

    //                 const size_t path2_start = path2.first.start_pos();
    //                 const size_t path1_end = path.first.end_pos();

    //                 if(path2.first.front().first == path.first.back().first){
    //                     long int distance_between_reads = path2_start - path1_end;
    //                     if(abs(distance_between_reads) < cfg::get().barcode_distance && distance_between_reads >= 0)
    //                         AddEdge(visited, path, path2, path_set, gp, first, barcode);
    //                 } else{
    //                     DistancesLengthsCallback<debruijn_graph::DeBruijnGraph> callback(gp.g);
    //                     VertexId startVertex = gp.g.EdgeEnd(path.first.back().first);
    //                     VertexId endVertex = gp.g.EdgeStart(path2.first.front().first);
    //                     ProcessPaths(gp.g, 0, cfg::get().barcode_distance, startVertex, endVertex, callback);
    //                     if(callback.distances().size()){
    //                         vector<size_t> distances = callback.distances();
    //                         size_t min_elem = distances[0];
    //                         for(size_t i = 1; i < distances.size(); ++i){
    //                             if(distances[i] < min_elem) min_elem = distances[i];
    //                         }
    //                         long int distance_between_reads = cfg::get().barcode_distance - path2_start - path1_end - min_elem;
    //                         if(distance_between_reads >= 0) AddEdge(visited, path, path2, path_set, gp, first, barcode);
    //                     }
    //                 }
    //             }
    //         }
    //         // for each component with more than one paired read output another barcode, such as -1, -2 etc
    //     }
    // }


    std::vector<VertexId> VerticesReachedFrom(VertexId& start_vertex, 
                        debruijn_graph::conj_graph_pack &gp, int edge_size) {
        // INFO("Vertices reached from distance: " << edge_size); // LM: Fix_2. This is the actual maximum path length of the Dijkstra graph.
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(start_vertex);
        // INFO("Reached vertices size - " << bounded_dijkstra.ReachedVertices()); // LM: Testing. TRACE() isn't showing up even with --debug.
        return bounded_dijkstra.ReachedVertices();
    }
    std::vector<VertexId> ConjugateVerticesReachedFrom(VertexId& start_vertex, 
                        debruijn_graph::conj_graph_pack &gp, int edge_size) {
        // INFO("vertices reached from distace: " << edge_size); // LM: Fix_2. This is the actual maximum path length of the backwards Dijkstra graph.
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(gp.g, edge_size);
        bounded_dijkstra.Run(start_vertex);
        // INFO("Reached vertices size - " << bounded_dijkstra.ReachedVertices()); // LM: Testing. TRACE() isn't showing up even with --debug.
        return bounded_dijkstra.ReachedVertices();
    }

    void clusterReads(debruijn_graph::conj_graph_pack &gp,
        std::vector<std::pair<MappingPath<EdgeId>, std::pair<std::pair<std::string, std::string>, std::string>>>& paths,
        std::unordered_map<std::string, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>>& path_set,
        std::string& barcode,
        std::queue<path_extend::BidirectionalPath*>& stability_queue){
        std::unordered_map<const MappingPath<EdgeId>*, path_extend::BidirectionalPath*> visited;
        // int path_ct = 0; // LM: Testing. For R2 unit test.
        for (auto const& path : paths) {
            /* if (path_ct > 2) break; // LM: Testing.
            ++path_ct; //LM: Testing. */
            bool first = true;
            VertexId startVertex = gp.g.EdgeEnd(path.first.back().first);
            // LM: Fix_2. The end vertex of the reverse complement of the first edge of the read.
            VertexId conjStartVertex = gp.g.EdgeStart(path.first.front().first); // This is the start vertex of the edge that the read starts on.
            std::vector<VertexId> reached_vertices;
            std::vector<VertexId> conjugate_reached_vertices;
            // INFO(path.first.start_pos() << " " << path.first.end_pos()); // LM: Testing. This should print the start and end positions of the first read.
            // INFO(path.first.mapping_at(0).mapped_range.start_pos << " " << path.first.mapping_at(path.first.size() - 1).mapped_range.end_pos); // LM: Testing. This should print the entire mapping path of the read. See method at mapping_path.hpp line 260.
            // Hopefully, end_pos_read_on_edge = range_mappings_[range_mappings_.size() - 1].mapped_range.end_pos)
            // LM: Fix_2. Takes the difference between the end of the edge that the read-end sits on and the position of the read-end on the edge, instead of the end and start of the read. Replace entire if-loop with this.
            int endDist = gp.g.length(path.first.edge_at(path.first.size() - 1)) - path.first.end_pos();
            int startDist = path.first.mapping_at(0).mapped_range.start_pos;
            int reducedEndDist = cfg::get().barcode_distance - endDist;
            int reducedStartDist = cfg::get().barcode_distance - startDist;
            // INFO("End_Edge_Length = " << gp.g.length(path.first.edge_at(path.first.size() - 1)));
            // INFO("End = " << endDist << " Start = " << startDist << " Reduced_FD = " << reducedEndDist << " Reduced_CD = " << reducedStartDist); // LM: Testing. See if the reductions to the barcode distance are accurate.
            // INFO("End_Vertex = " << startVertex << " Conj_Start_Vertex = "<< conjStartVertex); // LM: Testing. See if we're printing out different vertices.
            if (reducedEndDist > 0){
                reached_vertices = VerticesReachedFrom(startVertex, gp, reducedEndDist); // Find the list of vertices that can be reached with a path of MappingPath<EdgeId>.;
            }
            if (reducedStartDist > 0){
                conjugate_reached_vertices = ConjugateVerticesReachedFrom(conjStartVertex, gp, reducedStartDist);
            }

            std::sort(reached_vertices.begin(), reached_vertices.end());
            std::sort(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end());
            for(auto const& path2 : paths) {
                if(&path.first != &path2.first){ // TODO: See if this equivalency works.
                    // If Read hasn't been traversed before & corresponding bidirectional path has not been made

                    // int counter = 0; // LM: Testing. Seeing how many times each second read is being added.
                    // INFO("Distance = " << distance_between_reads << " path_end = " << path1_end << " path2_start = " << path2_start);
                    // LM: Testing. The number of edges covered (partially) by the read.
                    // INFO("path_size = " << path.first.size() << " path2_size = " << path2.first.size());
                    // Fix_3A: The only thing to check is if the first read ends on the same read as the second read starts. Got rid of the two for-loops checking every edge.
                    // Original Fix_3: Small efficiency fix. If the reads i) overlap or ii) are too distant anyways then there's no point in checking to see if they overlap.
                    // for(size_t i = 0;i < path.first.size(); ++i){
                        // for(size_t j = 0;j < path2.first.size(); ++j){
                            // LM: Testing. This should give us the edge ID of the ith and jth edges of reads 1 and 2 respectively.
                            // INFO(path.first[i].first << " " << path2.first[j].first);
                            // LM: Testing. This should give us the mapping ranges of the ith and jth edges of reads 1 and 2 respectively.
                            // INFO(path.first[i].second.mapped_range << " " << path2.first[j].second.mapped_range);
                            // Fix_3A: If the edge ID of the read 1 end and the read 2 start are the same...
                            if (path.first.edge_at(path.first.size() - 1) == path2.first.edge_at(0)) {
                            // if(path.first[i].first == path2.first[j].first){
                                // LM: Fix_5: Previously, this was in the wrong place as it worked whether or not the reads ENDED on the same edge (i.e.: it could work when they coincided).
                                long int distance_between_reads = path2.first.start_pos() - path.first.end_pos();
                                if(std::abs(distance_between_reads) < cfg::get().barcode_distance && distance_between_reads >= 0){ // LM: Orignally, this was std::abs(distance) for the second term, which would have rendered it always true.
                                    // INFO("Same edge end. This is connection number " << counter); // LM: Testing. To find out how many times each read-read connection was added to the list.
                                    // ++counter; // LM: Testing. Tracking the above.
                                    AddEdge(visited, path, path2, path_set, gp, barcode, first, stability_queue);
                                    goto added; // LM: Boolean done was not a clean way to record whether or not connectivity established.
                                } 
                            // Fix_3B: To allow for the case that the conjugate of read 2 lies on the same end-edge as read 1. 
                            } else if (path.first.edge_at(path.first.size() - 1) == gp.g.conjugate(path2.first.edge_at(path2.first.size() - 1))) {
                                int path2_conj_end = gp.g.length(path2.first.edge_at(path2.first.size() - 1)) - path2.first.end_pos();
                                long int distance_between_reads = path2_conj_end - path.first.end_pos();
                                if(std::abs(distance_between_reads) < cfg::get().barcode_distance && distance_between_reads >= 0){ // LM: Orignally, this was std::abs(distance) for the second term, which would have rendered it always true.
                                    // INFO("Same edge end. This is connection number " << counter); // LM: Testing. To find out how many times each read-read connection was added to the list.
                                    // ++counter; // LM: Testing. Tracking the above.
                                    AddEdge(visited, path, path2, path_set, gp, barcode, first, stability_queue);
                                    goto added; // LM: Boolean done was not a clean way to record whether or not connectivity established.
                                }
                            }
                        // }
                    // }
                    for(size_t i = 0; i < path2.first.size(); ++i){
                        VertexId endVertex = gp.g.EdgeEnd(path2.first[i].first);
                        if (std::binary_search(reached_vertices.begin(), reached_vertices.end(), endVertex) ||
                            std::binary_search(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end(), endVertex)){
                            // INFO("Different edge end."); // LM: Testing. To find out how many times each read-read connection was added to the list.
                            AddEdge(visited, path, path2, path_set, gp, barcode, first, stability_queue);
                        } else { // LM: Fix_1. Find the reverse complement of the second read, and check it against the set of (reverse complement) reached vertices.
                            endVertex = gp.g.conjugate(gp.g.EdgeEnd(path2.first[i].first)); // This is the start vertex of an edge that the read sits on. IOW, the end vertex of an edge that the read reverse complement sits on.
                            //  Fix_6: Previously, didn't have gp.g.conjugate() so it wasn't really the conjugate. 
                            if (std::binary_search(reached_vertices.begin(), reached_vertices.end(), endVertex) ||
                                std::binary_search(conjugate_reached_vertices.begin(), conjugate_reached_vertices.end(), endVertex)) {
                                // INFO("Different conjugate edge end"); // LM: Testing. To find out how many times each read-read connection was added to the list.
                                AddEdge(visited, path, path2, path_set, gp, barcode, first, stability_queue);
                            }
                        }
                    }
                    // DistancesLengthsCallback<debruijn_graph::DeBruijnGraph> callback(gp.g);
                    // VertexId startVertex = gp.g.EdgeEnd(path.first.back().first);
                    // VertexId endVertex = gp.g.EdgeStart(path2.first.front().first);
                    // ProcessPaths(gp.g, 0, cfg::get().barcode_distance, startVertex, endVertex, callback);
                    // if(callback.distances().size()){
                    //     vector<size_t> distances = callback.distances();
                    //     size_t min_elem = distances[0];
                    //     for(size_t i = 1; i < distances.size(); ++i){
                    //         if(distances[i] < min_elem) min_elem = distances[i];
                    //     }
                    // long int distance_between_reads = cfg::get().barcode_distance - path2_start - path1_end - min_elem;
                    // if(distance_between_reads >= 0) AddEdge(visited, path, path2, path_set, gp, first, barcode);
                    // }
                    
                }
                added:;
            }
        }
        for (auto const& path : paths) {
            if(!visited.count(&path.first)){
                std::vector<EdgeId> edges = path.first.simple_path();
                path_extend::BidirectionalPath* bidirectional_path = new path_extend::BidirectionalPath(gp.g);
                bidirectional_path->SetConjPath(new path_extend::BidirectionalPath(gp.g));
                for (auto e : edges) {
                    bidirectional_path->PushBack(e);
                    // Fix_4: TODO Removed GetConjPath() entirely because it's not reported at all, and there's really no point because the reverse conjugate connection technically doesn't exist.
                    bidirectional_path->GetConjPath()->PushBack(gp.g.conjugate(e));
                }
                // TODO -stop using pairs and store header(name), quality and sequence smarter
                bidirectional_path->name = path.second.first.first;
                bidirectional_path->quality_string = path.second.first.second;
                bidirectional_path->sequence_string = path.second.second;
                //                                          Pointer to path        vector of paths with edges to first(adjacencies)
                path_set[barcode][bidirectional_path] =  std::vector<path_extend::BidirectionalPath*>();
                stability_queue.push(bidirectional_path);
                visited[&path.first] = bidirectional_path;
                // first = false;
            }
        }
    }



    void processReads(debruijn_graph::conj_graph_pack &graph_pack, const lib_t& lib_10x) {
        auto mapper = MapperInstance(graph_pack);
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;
        std::vector<std::pair<MappingPath<EdgeId>, std::pair<std::pair<std::string, std::string>, std::string>>> paths;
        std::queue<path_extend::BidirectionalPath*> stability_queue;
        // INFO("barcode distance: " << cfg::get().barcode_distance);


        path_extend::FastqWriter writer2(graph_pack.g, make_shared<path_extend::DefaultContigNameGenerator>());
        std::string file_name = cfg::get().output_dir + std::to_string(cfg::get().barcode_distance) + "extracted.fastq";
        // INFO("Outputting updated reads with enhanced barcodes to " << file_name);
        io::OFastqReadStream os_(file_name);
        std::ofstream statistics_file;
        statistics_file.open(cfg::get().output_dir + std::to_string(cfg::get().barcode_distance) + "statistics.txt");


        //barcode --> BidirectionalPath
        std::unordered_map<std::string, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>> connected_components;
        int first_thousand = 0;
        std::string current_barcode = "";
        while(!stream->eof()) {
            *stream >> read;
            std::string barcode_string = GetTenXBarcodeFromRead(read);
            if(barcode_string != ""){
                // if(first_thousand >= 1000) break;
                if(barcode_string != current_barcode && !paths.empty()){ //  && first_thousand < 1000
                    first_thousand++;
                    if(current_barcode == ""){
                        current_barcode = GetTenXBarcodeFromRead(paths[0].second.first.first);
                    }
                    // LM: Testing. Skipping the first already-done barcodes so that it doesn't take forever.
                    if (first_thousand < 920225) {
                        INFO("Skipping barcodes " << current_barcode  << " number " << first_thousand);
                        paths.clear();
                        continue;
                    }
                    INFO(first_thousand << ": Processing barcode " << current_barcode << ": " << paths.size() << "(Number of reads in barcode)");
                    clusterReads(graph_pack, paths, connected_components, current_barcode, stability_queue);
                    stability_queue = writer2.OutputPaths(connected_components[current_barcode], current_barcode, os_, statistics_file, stability_queue);
                    int pewpew = connected_components.erase(current_barcode);
                    // LM: Testing. Stops unit test after first barcode with >10 reads, which will likely be separated into >5 enhanced barcodes.
                    /* if (paths.size() > 10) {
                        INFO("Should terminate at barcode " << current_barcode);
                        paths.clear();
                        break;
                    } */
                    paths.clear();
                }
                const auto &path1 = mapper->MapRead(read.first());
                const auto &path2 = mapper->MapRead(read.second());
                // LEFT OFF HERE
                if(path1.size()) {
                    paths.push_back(std::make_pair(path1, std::make_pair(std::make_pair(read.first().name(), read.first().GetQualityString()), read.first().GetSequenceString())));
                }
                if(path2.size()) {
                    paths.push_back(std::make_pair(path2, std::make_pair(std::make_pair(read.second().name(), read.second().GetQualityString()), read.second().GetSequenceString())));
                }
            }
            current_barcode = barcode_string;
        }
        /* if(first_thousand < 1000){
            first_thousand++;
            INFO(first_thousand << ": Processing barcode " << current_barcode << ": " << paths.size() << "(Number of reads in barcode)");
            clusterReads(graph_pack, paths, connected_components, current_barcode, stability_queue);
        } */
        stability_queue = writer2.OutputPaths(connected_components[current_barcode], current_barcode, os_, statistics_file, stability_queue);
        int pewpew = connected_components.erase(current_barcode);
        paths.clear();
    }
}

#endif
