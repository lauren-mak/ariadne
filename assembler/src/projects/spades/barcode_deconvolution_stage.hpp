#ifndef BARCODE_DECONVOLUTION_STAGE_HPP
#define BARCODE_DECONVOLUTION_STAGE_HPP

#include "pipeline/stage.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include <unordered_map>
#include <string>


namespace debruijn_graph {

    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> lib_t;
    typedef std::vector<io::SequencingLibrary<debruijn_graph::config::LibraryData>> lib_vector_t;

    class SimplePathExtractor;
    class LongReadsCreator;
    std::string GetTenXBarcodeFromRead(const io::PairedRead &read);
    void processReads(debruijn_graph::conj_graph_pack &graph_pack, const lib_t& lib_10x);

    class BarcodeDeconvolutionStage : public spades::AssemblyStage {
    public:

        BarcodeDeconvolutionStage() : AssemblyStage("10X dbg", "10X_dbg") {}

      	void run(conj_graph_pack &gp, const char*){
            gp.EnsureIndex();
            if(!gp.kmer_mapper.IsAttached()) gp.kmer_mapper.Attach();
            INFO("Barcode Deconvolution Started");
            config::dataset& dataset_info = cfg::get_writable().ds;
            for (auto& lib : dataset_info.reads) {
                INFO("START PROCESSING READS");
                processReads(gp, lib);
            }
            INFO("Barcode Deconvolution Ended");
        }

    };


    void extendBackward(debruijn_graph::conj_graph_pack &gp, const std::set<EdgeId> &edge_set, std::deque<EdgeId> &linear_path, EdgeId e, std::set<EdgeId> &used_edges);


    std::vector<EdgeId> getForwardIntersect(debruijn_graph::conj_graph_pack &gp, EdgeId e, const std::set<EdgeId> &edge_set, std::set<EdgeId>& used_edges) {
        auto edges = gp.g.IncidentEdges(gp.g.EdgeEnd(e));
        std::vector<EdgeId> filtered;
        for (auto temp_e : edges) {
            if (temp_e != e && edge_set.count(temp_e) && !used_edges.count(temp_e)) {
                filtered.push_back(temp_e);
            }
        }
        return filtered;
    }

    std::vector<EdgeId> getReverseIntersect(debruijn_graph::conj_graph_pack &gp, EdgeId e, const std::set<EdgeId> &edge_set, std::set<EdgeId>& used_edges) {
        auto edges = gp.g.IncidentEdges(gp.g.EdgeStart(e));
        std::vector<EdgeId> filtered;
        for (auto temp_e : edges) {
            if (temp_e != e && edge_set.count(temp_e) && !used_edges.count(temp_e)) {
                filtered.push_back(temp_e);
            }
        }
        return filtered;
    }

    std::vector<EdgeId> getForwardIntersect(debruijn_graph::conj_graph_pack &gp, EdgeId e, const GraphComponent<Graph> &comp) {

        auto edges = gp.g.IncidentEdges(gp.g.EdgeEnd(e));
        DEBUG("Edges " << edges);
        std::vector<EdgeId> filtered;
        for (auto temp_e : edges) {
            if (temp_e != e && comp.edges().count(temp_e)) {
                filtered.push_back(temp_e);
            }
        }
        DEBUG("Filtered " << filtered);
        return filtered;
    }

    std::vector<EdgeId> getReverseIntersect(debruijn_graph::conj_graph_pack &gp, EdgeId e, const GraphComponent<Graph> &comp) {
        auto edges = gp.g.IncidentEdges(gp.g.EdgeStart(e));
        std::vector<EdgeId> filtered;
        for (auto temp_e : edges) {
            if (temp_e != e && comp.edges().count(temp_e)) {
                filtered.push_back(temp_e);
            }
        }
        return filtered;
    }

    void extendForward(debruijn_graph::conj_graph_pack &gp, const GraphComponent<Graph> &comp,
                       std::deque<EdgeId> &linear_path, std::set<EdgeId> &used_edges) {
        EdgeId e = linear_path.back();
        if (gp.g.IncomingEdgeCount(gp.g.EdgeEnd(e)) != 1 ||
                gp.g.OutgoingEdgeCount(gp.g.EdgeEnd(e)) != 1 ) {
            return;
        }
        auto extensions = getForwardIntersect(gp, e, comp);
        VERIFY(extensions.size() == 1);
        EdgeId next_edge = extensions[0];
        if (used_edges.count(next_edge)) {
            return;
        }
        linear_path.push_back(next_edge);
        used_edges.insert(next_edge);
        used_edges.insert(gp.g.conjugate(next_edge));
        extendForward(gp, comp, linear_path, used_edges);
    }

    void extendBackward(debruijn_graph::conj_graph_pack &gp, const GraphComponent<Graph> &comp,
                        std::deque<EdgeId> &linear_path, std::set<EdgeId> &used_edges) {
        EdgeId e = linear_path.front();
        if (gp.g.IncomingEdgeCount(gp.g.EdgeStart(e)) != 1 ||
                gp.g.OutgoingEdgeCount(gp.g.EdgeStart(e)) != 1 ) {
            return;
        }
        auto extensions = getReverseIntersect(gp, e, comp);
        VERIFY(extensions.size() == 1);
        EdgeId prev_edge = extensions[0];
        if (used_edges.count(prev_edge)) {
            return;
        }
        linear_path.push_front(prev_edge);
        used_edges.insert(prev_edge);
        used_edges.insert(gp.g.conjugate(prev_edge));
        extendBackward(gp, comp, linear_path, used_edges);
    }

    void extendForward(debruijn_graph::conj_graph_pack &gp, const std::set<EdgeId> &edge_set,
                       std::deque<EdgeId> &linear_path, EdgeId e, std::set<EdgeId>& used_edges) {
        auto extensions = getForwardIntersect(gp, e, edge_set, used_edges);
        for (auto next_edge : extensions) {
            linear_path.push_back(next_edge);
            used_edges.insert(e);
            used_edges.insert(gp.g.conjugate(e));
            extendForward(gp, edge_set, linear_path, next_edge, used_edges);
            extendBackward(gp, edge_set, linear_path, next_edge, used_edges);
        }
    }


    void extendBackward(debruijn_graph::conj_graph_pack &gp, const std::set<EdgeId> &edge_set,
                        std::deque<EdgeId> &linear_path, EdgeId e, std::set<EdgeId>& used_edges) {
        auto extensions = getReverseIntersect(gp, e, edge_set, used_edges);
        for (auto prev_edge : extensions) {
            linear_path.push_front(prev_edge);
            used_edges.insert(prev_edge);
            used_edges.insert(gp.g.conjugate(prev_edge));
            extendForward(gp, edge_set, linear_path, prev_edge, used_edges);
            extendBackward(gp, edge_set, linear_path, prev_edge, used_edges);
        }
    }

    bool IsSimplePath(const GraphComponent<Graph> &comp) {
        size_t one_degree = 0;
        for (auto v : comp.vertices()) {
            auto edges = comp.IncidentEdges(v);
            VERIFY(edges.size() > 0);
            if (edges.size() >= 3) {
                return false;
            }
            if (edges.size() == 1) {
                one_degree++;
            }
        }
        return one_degree == 4;
    }

    //Here we suppose that this is single connected component
    bool IsSimpleCycle(const GraphComponent<Graph> &comp) {
        return comp.vertices().size() == 2 && comp.edges().size() == 2;
    }

    size_t SplitComponentsOnSimplePaths(debruijn_graph::conj_graph_pack &gp, const GraphComponent<Graph> &comp, path_extend::PathContainer &temp_set) {
        size_t result = 0;
        //check self-conjugate
        set<EdgeId> used_edges;
        for (auto e : comp.edges()) {
            if (gp.g.conjugate(e) == e) {
                used_edges.insert(e);
            }
            if (gp.g.EdgeEnd(e) == gp.g.EdgeStart(e)) {
                used_edges.insert(e);
            }
        }


        DEBUG(comp.edges());
        for (auto e : comp.edges()) {
            if (!used_edges.count(e)) {
                DEBUG("Current edge " << e);
                used_edges.insert(e);
                used_edges.insert(gp.g.conjugate(e));
                std::deque<EdgeId> linear_path;
                linear_path.push_back(e);
                extendForward(gp, comp, linear_path, used_edges);
                extendBackward(gp, comp, linear_path, used_edges);
                path_extend::BidirectionalPath *path = new path_extend::BidirectionalPath(gp.g);
                auto conj = new path_extend::BidirectionalPath(gp.g);
                temp_set.AddPair(path, conj);
                for (auto e : linear_path) {
                    path->PushBack(e);
                }
                DEBUG("Add path");
                path->PrintDEBUG();
                result++;
            }
        }
        return result;
    }

    void analogOfExtractLongReads(debruijn_graph::conj_graph_pack &gp, std::vector<MappingPath<EdgeId>>& paths, path_extend::PathContainer& path_set, std::string& barcode) {
        std::set<EdgeId> edge_set;
        for (auto const& path : paths) {
            std::vector<EdgeId> edges = path.simple_path();
            for (auto e : edges) {
                edge_set.insert(e);
                edge_set.insert(gp.g.conjugate(e));
            }
        }
        path_extend::PathContainer temp_set;



        std::set<EdgeId> used_edges;
        for (auto e : edge_set) {
            if (!used_edges.count(e)) {

                std::deque<EdgeId> linear_path;
                linear_path.push_back(e);
                used_edges.insert(e);
                used_edges.insert(gp.g.conjugate(e));
                extendForward(gp, edge_set, linear_path, e, used_edges);
                extendBackward(gp, edge_set, linear_path, e, used_edges);
                auto component = GraphComponent<Graph>::FromEdges(gp.g, linear_path, true);
                // component.ClipTips();

                if (IsSimplePath(component)) {
                    DEBUG("Component is a simple path");
                    DEBUG(component.edges());
                    path_extend::BidirectionalPath *path = new path_extend::BidirectionalPath(gp.g);
                    path->SetConjPath(new path_extend::BidirectionalPath(gp.g));
                    for (auto e : linear_path) {
                        path->PushBack(e);
                    }
                    temp_set.AddPair(path, path->GetConjPath());
                } else if (IsSimpleCycle(component)) {
                    DEBUG("Component is a simple cycle");
                    DEBUG(component.edges());
                } else {
                    DEBUG("Component is not a simple path");
                    DEBUG(component.edges());
                    size_t result = SplitComponentsOnSimplePaths(gp, component, temp_set);
                    DEBUG("Component is split on " << result << " paths");
                }
            }
        }


        for (auto path_pair : temp_set) {
            // INFO("DO I MAKE IT HERE " << barcode);
            path_pair.first->barcode = barcode;
            path_pair.second->barcode = barcode;

            path_set.AddPair(path_pair.first, path_pair.second);
        }
        temp_set.clear();
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
            TRACE(barcode);
            return barcode;
        }
        return "";
    }


    void processReads(debruijn_graph::conj_graph_pack &graph_pack, const lib_t& lib_10x) {
        INFO("CHECK PROCESS READS STARTS");
        auto mapper = MapperInstance(graph_pack);
        auto stream = io::paired_easy_reader(lib_10x, false, false);
        io::PairedRead read;
        std::vector<MappingPath<EdgeId>> paths;
        std::string current_barcode = "";
        // path_extend::PathContainer& long_reads = graph_pack.mapped_paths[lib_10x.data().lib_index];
        // ^^ Originally this, no member mapped_paths in graph_pack, maybe contig_paths
        path_extend::PathContainer long_reads;
        // LongReadsCreator extractor(graph_pack);
        while(!stream->eof()) {
            *stream >> read;
            string barcode_string = GetTenXBarcodeFromRead(read);
            if(barcode_string != ""){
                if(barcode_string != current_barcode && !paths.empty()) {
                  // long_reads = // vector of pairs(pointer to path(deque of edgeids),
                  // pointer to path for conjugate graph(deque of edgeids))
                  analogOfExtractLongReads(graph_pack, paths, long_reads, current_barcode);
                  paths.clear();
                }
            }
            const auto &path1 = mapper->MapRead(read.first());
            const auto &path2 = mapper->MapRead(read.second());
            paths.push_back(path1);
            paths.push_back(path2);

            current_barcode = barcode_string;
        }
        analogOfExtractLongReads(graph_pack, paths, long_reads, current_barcode);
        paths.clear();
        bool distances_found = false;

        // struct hash_pairs{
        //     std::size_t operator() (const std::pair<std::string, std::string>& p) const {
        //         return boost::hash_value(p);
        //     }
        // };
        // std::unordered_map<std::string, path_extend::PathContainer> associated_barcodes;
        std::unordered_map<std::string, std::vector<size_t>> associated_barcodes;
        

        for(auto pair1 : long_reads) {
            for(auto pair2 : long_reads){
                

                // if(&pair1.first != &pair2.first && pair1.first->barcode == pair2.first->barcode) {
                //     VertexId startVertex = graph_pack.g.EdgeEnd(pair1.first->Back());
                //     VertexId endVertex = graph_pack.g.EdgeStart(pair2.first->Front());
                //     // I have the vertices, just find distance. ProcessPath gives you the path
                //     DistancesLengthsCallback<debruijn_graph::DeBruijnGraph> callback(graph_pack.g);
                //     ProcessPaths(graph_pack.g, 0, 25000, startVertex, endVertex, callback);
                //     std::vector<size_t> all_distances = callback.distances();
                //     if(all_distances.size() > 0) {
                //         distances_found = true;
                //         if(associated_barcodes.find(pair1.first->barcode) != associated_barcodes.end()){
                //             associated_barcodes[pair1.first->barcode].AddPair(pair1.first, pair1.second);
                //             associated_barcodes[pair1.first->barcode].AddPair(pair2.first, pair2.second);
                //         } else{
                //             associated_barcodes[pair1.first->barcode];
                //             associated_barcodes[pair1.first->barcode].AddPair(pair1.first, pair1.second);
                //             associated_barcodes[pair1.first->barcode].AddPair(pair2.first, pair2.second);
                //         }
                //     }



                //     for(auto howFar : all_distances){
                //         INFO("distance from " << pair1.first->barcode << " " << pair1.first->GetId() << " to " << pair2.first->barcode << " " << pair2.first->GetId() << ": " << howFar);
                //     }
                // }
                if(&pair1.first != &pair2.first && pair1.first->barcode == pair2.first->barcode) {
                    VertexId startVertex = graph_pack.g.EdgeEnd(pair1.first->Back());
                    VertexId endVertex = graph_pack.g.EdgeStart(pair2.first->Front());
                    // I have the vertices, just find distance. ProcessPath gives you the path
                    DistancesLengthsCallback<debruijn_graph::DeBruijnGraph> callback(graph_pack.g);
                    ProcessPaths(graph_pack.g, 0, 25000, startVertex, endVertex, callback);
                    std::vector<size_t> all_distances = callback.distances();
                    if(all_distances.size() > 0) {
                        distances_found = true;
                        if(associated_barcodes.find(pair1.first->barcode) != associated_barcodes.end()){
                            associated_barcodes[pair1.first->barcode].push_back(pair1.first->GetId());
                            associated_barcodes[pair1.first->barcode].push_back(pair2.first->GetId());
                        } else{
                            associated_barcodes[pair1.first->barcode];
                            associated_barcodes[pair1.first->barcode].push_back(pair1.first->GetId());
                            associated_barcodes[pair1.first->barcode].push_back(pair2.first->GetId());
                        }
                    }



                    for(auto howFar : all_distances){
                        INFO("distance from " << pair1.first->barcode << " " << pair1.first->GetId() << " to " << pair2.first->barcode << " " << pair2.first->GetId() << ": " << howFar);
                    }
                }
            }
        }
        for(auto it = associated_barcodes.begin(); it != associated_barcodes.end(); ++it) {
            INFO(it->first << ": " << it->second.size());
            path_extend::ContigWriter writer(graph_pack.g, make_shared<path_extend::DefaultContigNameGenerator>());
            INFO("Outputting individual barcode clusters for: " << it->first << " to " << cfg::get().output_dir << "extracted.fasta");
            // writer.OutputPaths(it->second, cfg::get().output_dir + "barcodes/" + it->first + "extracted.fasta");
        }
        // path_extend::ContigWriter writer(graph_pack.g, make_shared<path_extend::DefaultContigNameGenerator>());
        // INFO("Outputting updated reads with barcode to " << cfg::get().output_dir << "extracted.fasta");
        // writer.OutputPaths(long_reads, cfg::get().output_dir + "extracted.fasta");
        // boost::disjoint_sets
    }

}

#endif