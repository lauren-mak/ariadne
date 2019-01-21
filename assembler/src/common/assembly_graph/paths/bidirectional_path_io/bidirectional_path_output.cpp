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

// void depthFirstSearch(BidirectionalPath*& first, std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>& paths, std::unordered_set<BidirectionalPath*>& visited, io::OFastqReadStream& os_, size_t& counter, size_t& statistics){
//     for(size_t i = 0; i < paths[first].size(); ++i){
//         if(!visited.count(paths[first][i])){
//             std::string name = paths[first][i]->name;
//             std::string sequence_string = paths[first][i]->sequence_string;
//             std::string quality_string = paths[first][i]->quality_string;
//             io::SingleRead left(name + " cluster#: " + std::to_string(counter), sequence_string, quality_string);
//             os_ << left;
//             ++statistics;
//             visited.insert(paths[first][i]);
//             depthFirstSearch(paths[first][i], paths, visited, os_, counter, statistics);
//         }
//     }
// }

std::string GetUpdatedReadName(std::string& read, size_t& counter) {
    std::string delimeter = "-1";
    size_t start_pos = read.find(delimeter);
    size_t delimeter_size = delimeter.length();
    if (start_pos != string::npos) {
        std::string read_name = 
            read.substr(0, start_pos) + "-" + std::to_string(counter) + read.substr(start_pos+delimeter_size);
        TRACE(read_name);
        return read_name;
    }
    return "";
}

void path_extend::FastqWriter::OutputPaths(std::unordered_map<path_extend::BidirectionalPath*, std::vector<path_extend::BidirectionalPath*>>& paths, 
    std::string& barcode, 
    io::OFastqReadStream& os_, 
    std::ofstream& statistics_file) {
    
    std::queue<BidirectionalPath*> cluster_queue;

    std::unordered_set<BidirectionalPath*> visited;
    // INFO("Number of nodes in " << it->first << ": " << paths[it->first].size());

    size_t counter = 1;
    for(auto iter = paths.begin(); iter != paths.end(); ++iter){
        if(visited.count(iter->first)) continue;
        size_t statistics = 0;
        
        std::string name = GetUpdatedReadName(iter->first->name, counter);
        std::string sequence_string = iter->first->sequence_string;
        std::string quality_string = iter->first->quality_string;

        // read from which connected component graph is built
        //running a BFS

        if(paths[iter->first].size() > 1 && !visited.count(iter->first)){

            cluster_queue.push(iter->first);
            visited.insert(iter->first);
            while(!cluster_queue.empty()){
                // x:dequeue()
                BidirectionalPath* pointer = cluster_queue.front();
                
                // Mark x in T
                std::string name = GetUpdatedReadName(iter->first->name, counter);
                std::string sequence_string = pointer->sequence_string;
                std::string quality_string = pointer->quality_string;
                io::SingleRead left(name, sequence_string, quality_string);
                os_ << left;
                ++statistics;
                
                // For all unmarked neighbors of X
                for (size_t i = 0; i < paths[pointer].size(); ++i){
                    if(!visited.count(paths[pointer][i])){
                        // Mark and put in queue
                        cluster_queue.push( paths[pointer][i] );
                        visited.insert( paths[pointer][i] );
                    }
                }
                cluster_queue.pop();
            }
            statistics_file << barcode << "-" << counter << ", number reads in cluster: " << statistics << endl;

            // // all reads in singular connected component
            // for(size_t i = 0; i < paths[iter->first].size(); ++i){

            //     if(!visited.count(paths[iter->first][i])){
            //         std::string other_name = paths[iter->first][i]->name;
                    
            //         std::string other_sequence_string = paths[iter->first][i]->sequence_string;
            //         std::string other_quality_string = paths[iter->first][i]->quality_string;
            //         io::SingleRead other(other_name + " cluster#: " + std::to_string(counter), other_sequence_string, other_quality_string);
            //         os_ << other;
            //         ++statistics;
            //         visited.insert(paths[iter->first][i]);
            //         depthFirstSearch(paths[iter->first][i], paths, visited, os_, counter, statistics);
            //         ++counter;
            //     }
            // }
            
        } 
        else if (!visited.count(iter->first)){
            visited.insert(iter->first);
            io::SingleRead left(name, sequence_string, quality_string);
            statistics_file << barcode << "-" << counter <<  ": singleton cluster" << endl;
            os_ << left;
            ++statistics;
        }
        counter++;
    }
}

}