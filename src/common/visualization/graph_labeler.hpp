//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/stl_utils.hpp"
#include <common/barcode_index/barcode_index.hpp>
#include <common/barcode_index/barcode_info_extractor.hpp>
#include "utils/standard_base.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"

namespace visualization {

namespace graph_labeler {

/**
* (Interface)
* Provides string labels for vertices and edges of some graph.
* Used with GraphPrinter to visualize graphs.
*/
template<class Graph>
class GraphLabeler {
public:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    virtual ~GraphLabeler() {
    }

    virtual string label(VertexId v) const = 0;

    virtual string label(EdgeId e) const = 0;

};

//template<class Graph>
//class MapGraphLabeler {
//    typedef typename Graph::EdgeId EdgeId;
//    typedef typename Graph::VertexId VertexId;
//    map<EdgeId, string> edge_map_;
//    map<VertexId, string> vertex_map_;
//
//public:
//
//    string label(VertexId v) const {
//        auto it = vertex_map_.find(v);
//        if (it == vertex_map_.end())
//            return "";
//        else
//            return it->second;
//    }
//
//    string label(EdgeId e) const {
//        auto it = edge_map_.find(e);
//        if (it == edge_map_.end())
//            return "";
//        else
//            return it->second;
//    }
//
//};

template<class Graph>
class AbstractGraphLabeler : public GraphLabeler<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &g_;
protected:
    AbstractGraphLabeler(const Graph &g) : g_(g) {

    }

    const Graph &graph() const {
        return g_;
    }

public:
    /*virtual*/ std::string label(VertexId /*v*/) const {
        return "";
    }

    /*virtual*/ std::string label(EdgeId /*e*/) const {
        return "";
    }

};

/**
* Trivial implementation of GraphLabeler.
* All labels are "".
*/
template<class Graph>
class EmptyGraphLabeler : public GraphLabeler<Graph> {
    typedef GraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    EmptyGraphLabeler() {}

    std::string label(VertexId /*v*/) const {
        return "";
    }

    std::string label(EdgeId /*e*/) const {
        return "";
    }
};

/**
* Implementation of GraphLabeler for Graphs that have methods
* str(VertexId) and str(EdgeId), such as AbstractGraph.
*/
template<class Graph>
class StrGraphLabeler : public AbstractGraphLabeler<Graph> {
    typedef AbstractGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    StrGraphLabeler(const Graph &g) : base(g) {}

    /*virtual*/ std::string label(VertexId v) const {
        return this->graph().str(v);
    }

    /*virtual*/ std::string label(EdgeId e) const {
        return this->graph().str(e);
    }

    /*virtual*/ ~StrGraphLabeler() {

    }
};

template<class Graph>
shared_ptr<GraphLabeler<Graph>> StrGraphLabelerInstance(const Graph &g) {
    return make_shared<StrGraphLabeler<Graph>>(g);
}

template<class Graph>
class LengthIdGraphLabeler : public StrGraphLabeler<Graph> {
    typedef StrGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    LengthIdGraphLabeler(const Graph &g) : base(g) {}

    /*virtual*/ std::string label(EdgeId e) const {
        std::stringstream ss;
        ss << this->graph().length(e) << " (id: " << this->graph().int_id(e) << ")";
        return ss.str();
    }

};

template<class Graph>
class LengthGraphLabeler : public StrGraphLabeler<Graph> {
    typedef StrGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    LengthGraphLabeler(const Graph &g) : base(g) {}

    /*virtual*/ std::string label(EdgeId e) const {
        return std::to_string(this->graph().length(e));
    }

};

template<class Graph>
class CoverageGraphLabeler : public AbstractGraphLabeler<Graph> {
    typedef AbstractGraphLabeler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    CoverageGraphLabeler(const Graph &g) : base(g) {}

    std::string label(EdgeId e) const {
        double coverage = this->graph().coverage(e);
        return " {Cov:" + std::to_string(coverage) + "}";
    }
};

template<class Graph>
class CompositeLabeler : public GraphLabeler<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    vector<GraphLabeler<Graph> *> list_;

    template<typename ElementId>
    string ConstructLabel(ElementId id) const {
        vector<string> to_print;
        for (size_t i = 0; i < list_.size(); i++) {
            string next = list_[i]->label(id);
            if (next.size() != 0) {
                to_print.push_back(next);
            }
        }
        string result = "";
        for (size_t i = 0; i < to_print.size(); i++) {
            result += to_print[i];
            if (i + 1 < to_print.size())
                result += "\\n";
        }
        return result;
    }

public:
    CompositeLabeler() {
    }

    CompositeLabeler(GraphLabeler<Graph> &labeler1, GraphLabeler<Graph> &labeler2,
                     GraphLabeler<Graph> &labeler3,
                     GraphLabeler<Graph> &labeler4) {
        AddLabeler(labeler1);
        AddLabeler(labeler2);
        AddLabeler(labeler3);
        AddLabeler(labeler4);
    }

    CompositeLabeler(GraphLabeler<Graph> &labeler1, GraphLabeler<Graph> &labeler2,
                     GraphLabeler<Graph> &labeler3) {
        AddLabeler(labeler1);
        AddLabeler(labeler2);
        AddLabeler(labeler3);
    }

    CompositeLabeler(GraphLabeler<Graph> &labeler1, GraphLabeler<Graph> &labeler2) {
        AddLabeler(labeler1);
        AddLabeler(labeler2);
    }

    virtual ~CompositeLabeler() {
    }

    void AddLabeler(GraphLabeler<Graph> &labeler) {
        list_.push_back(&labeler);
    }

    virtual string label(VertexId vertexId) const {
        return ConstructLabel<VertexId>(vertexId);
    }

    virtual string label(EdgeId edgeId) const {
        return ConstructLabel<EdgeId>(edgeId);
    }
};

template<class Graph>
class EdgePosGraphLabeler : public AbstractGraphLabeler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
public:
    const omnigraph::EdgesPositionHandler<Graph> &edge_pos_;

    EdgePosGraphLabeler(const Graph &g, const omnigraph::EdgesPositionHandler<Graph> &edge_pos) :
            AbstractGraphLabeler<Graph>(g), edge_pos_(edge_pos) {
    }

    virtual std::string label(EdgeId edgeId) const {
        return "Positions: " + edge_pos_.str(edgeId);
    }

    virtual ~EdgePosGraphLabeler() {
//        TRACE("~EdgePosGraphLabeler");
    }

private:
    DECL_LOGGER("EdgePosGraphLabeler")
};

template<class Graph>
class DefaultLabeler : public GraphLabeler<Graph> {
private:
    const Graph &g_;
    const omnigraph::EdgesPositionHandler<Graph> &edges_positions_;
protected:
    typedef GraphLabeler<Graph> super;
    typedef typename super::EdgeId EdgeId;
    typedef typename super::VertexId VertexId;
public:

    DefaultLabeler(const Graph &g, const omnigraph::EdgesPositionHandler<Graph> &position_handler) :
            g_(g), edges_positions_(position_handler) {
    }

    virtual std::string label(VertexId vertexId) const {
        return std::to_string(vertexId.int_id());
    }

    virtual std::string label(EdgeId edgeId) const {
        std::string ret_label;
        ret_label += "Id " + g_.str(edgeId) + "\\n";
        ret_label += "Positions:\\n" + edges_positions_.str(edgeId);
        size_t len = g_.length(edgeId);
        double cov = g_.coverage(edgeId);
        ret_label += "Len(cov): " + std::to_string(len) + "(" + std::to_string(cov) + ")";
        return ret_label;
    }

    virtual ~DefaultLabeler() {
    }
};

template<class Graph>
class ReadCloudLabeler : public visualization::graph_labeler::StrGraphLabeler<Graph> {
    typedef visualization::graph_labeler::StrGraphLabeler <Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef barcode_index::FrameBarcodeIndexInfoExtractor extractor_t;
    typedef std::set <EdgeId> edge_set_t;
private:
    shared_ptr<extractor_t> barcode_extractor_ptr_;
    const size_t barcodes_on_string_;
    const size_t max_strings_;
    const size_t too_many_barcodes_;
public:
    ReadCloudLabeler(const Graph &g) :
            base(g), barcodes_on_string_(6), max_strings_(3), too_many_barcodes_(50) {

    }

    //fixme refactor barcode index to avoid this
    void UpdateExtractor(shared_ptr<barcode_index::AbstractBarcodeIndex> index_ptr, const Graph& g) {
        barcode_extractor_ptr_ = make_shared<extractor_t>(index_ptr, g);
    }

    virtual ~ReadCloudLabeler() {
    }

    std::string label(VertexId) const {
        return "";
    }

    std::string label(EdgeId e) const {
        std::string ret_label;
        size_t number_of_barcodes = barcode_extractor_ptr_->GetNumberOfBarcodes(e);
        vector <barcode_index::BarcodeId> barcodes = barcode_extractor_ptr_->GetBarcodes(e);
        size_t max_barcodes = barcodes_on_string_ * max_strings_;
        VERIFY(barcodes.size() == number_of_barcodes);
        ret_label += std::to_string(number_of_barcodes) + " barcodes.";
        if (number_of_barcodes > too_many_barcodes_) {
            return ret_label;
        }
        for (size_t index = 0; index < barcodes.size(); ++index) {
            if (index >= max_barcodes) {
                ret_label += "\\n and " + std::to_string(barcodes.size() - max_barcodes + 1) + " more.";
                break;
            }
            else {
                if (index % barcodes_on_string_ == 0) ret_label += "\\n";
                else ret_label += ", ";
                ret_label += std::to_string(barcodes[index]);
            }
        }
        return ret_label;
    }
};
}
}

