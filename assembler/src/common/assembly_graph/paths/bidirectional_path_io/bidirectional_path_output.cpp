//
// Created by andrey on 20.01.17.
//

#include "bidirectional_path_output.hpp"

namespace path_extend {

void path_extend::ContigWriter::OutputPaths(const PathContainer &paths, const vector<PathsWriterT>& writers) const {
    ScaffoldStorage storage;

    ScaffoldSequenceMaker scaffold_maker(g_);
    for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        BidirectionalPath* path = iter.get();
        if (path->Length() <= 0)
            continue;
        string path_string = scaffold_maker.MakeSequence(*path);
        if (path_string.length() >= g_.k()) {
            storage.emplace_back(path_string, path);
        }
    }

    //sorting by length and coverage
    std::sort(storage.begin(), storage.end(), [] (const ScaffoldInfo &a, const ScaffoldInfo &b) {
        if (a.length() == b.length())
            return math::gr(a.coverage(), b.coverage());
        return a.length() > b.length();
    });

    name_generator_->Preprocess(paths);
    for (size_t i = 0; i < storage.size(); ++i) {
        storage[i].name = name_generator_->MakeContigName(i+1, storage[i]);
    }

    for (auto& writer : writers) {
        writer(storage);
    }
    DEBUG("Contigs written");
}

void depthFirstSearch(BidirectionalPath*& first, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>& paths, std::unordered_set<BidirectionalPath*>& visited, io::OFastqReadStream& os_, size_t counter, const std::string& barcode, size_t& statistics){
    for(size_t i = 0; i < paths[first].size(); ++i){
        if(!visited.count(paths[first][i])){
            std::string sequence_string = paths[first][i]->sequence_string;
            std::string quality_string = paths[first][i]->quality_string;
            io::SingleRead left(barcode + "-" + std::to_string(counter), sequence_string, quality_string);
            os_ << left;
            ++statistics;
            visited.insert(paths[first][i]);
            depthFirstSearch(paths[first][i], paths, visited, os_, counter, barcode, statistics);
        }
    }
}

void path_extend::FastqWriter::OutputPaths(std::unordered_map<std::string, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>>& paths, std::string& file_){
    io::OFastqReadStream os_(file_);
    ofstream statistics_file;
    statistics_file.open(cfg::get().output_dir + "statistics.txt");
    for(auto barcode = paths.begin(); barcode != paths.end(); ++barcode){
        size_t counter = 1;
        std::unordered_set<BidirectionalPath*> visited;
        // INFO("Number of nodes in " << it->first << ": " << paths[it->first].size());
        for(auto iter = paths[barcode->first].begin(); iter != paths[barcode->first].end(); ++iter){
            size_t statistics = 0;
            std::string sequence_string = iter->first->sequence_string;
            std::string quality_string = iter->first->quality_string;
            if(paths[barcode->first][iter->first].size() > 1){
                if(!visited.count(iter->first)){
                    io::SingleRead left(barcode->first + "-" + std::to_string(counter), sequence_string, quality_string);
                    os_ << left;
                    ++statistics;
                    visited.insert(iter->first);
                }
                for(size_t i = 0; i < paths[barcode->first][iter->first].size(); ++i){

                    if(!visited.count(paths[barcode->first][iter->first][i])){
                        std::string other_sequence_string = paths[barcode->first][iter->first][i]->sequence_string;
                        std::string other_quality_string = paths[barcode->first][iter->first][i]->quality_string;
                        io::SingleRead other(barcode->first + "-" + std::to_string(counter), other_sequence_string, other_quality_string);
                        os_ << other;
                        ++statistics;
                        visited.insert(paths[barcode->first][iter->first][i]);
                        depthFirstSearch(paths[barcode->first][iter->first][i], paths[barcode->first], visited, os_, counter, barcode->first, statistics);
                        statistics_file << barcode->first << "-" << counter << ": " << statistics;
                        ++counter;
                    }
                }
                
            } else {
                io::SingleRead left(barcode->first, sequence_string, quality_string);
                os_ << left;
                ++statistics;
            }
            //barcodes, sizes of connected
        }
    }
}

}