#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>

using std::cout;
using std::endl;
using std::string;

const int MINOVERLAP = 10;

using Graph = boost::adjacency_list<
        boost::vecS,                                    // Edge container (std::vector for adjacency list)
        boost::vecS,                                    // Vertex container (std::vector for vertex list)
        boost::directedS,                               // Directed graph
        boost::property<boost::vertex_index_t, string>, // Vertex properties: integer id
        boost::property<boost::edge_weight_t, int>   // Edge properties integer weights (levenstien distance)
    >;


// Define a structure to hold each record
struct Record {
    std::string id;
    std::string pac;
    int start;
    int end;
    std::string sequence;
};

void print_record(Record r, bool print_seq = false){
    cout <<"id =" << r.id << " | start = " << r.start << " | end =" << r.end;
    if( print_seq ){
        cout << " | seq =  " << r.sequence;
    }
    cout << endl;
}

void print_interval(std::vector<int> I){
    cout<< "["<<I[0] << ", "<<I[1] <<"]" << endl;
}


int levenshteinDistance(const std::string& str1, const std::string& str2) {
    size_t len1 = str1.size();
    size_t len2 = str2.size();

    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1));

    for (size_t i = 0; i <= len1; ++i) {
        dp[i][0] = i; 
    }
    for (size_t j = 0; j <= len2; ++j) {
        dp[0][j] = j; // Cost of inserting characters into str1
    }

    // Fill the DP table
    for (size_t i = 1; i <= len1; ++i) {
        for (size_t j = 1; j <= len2; ++j) {
            // If characters are the same, no cost
            int cost = (str1[i - 1] == str2[j - 1]) ? 0 : 1;

            // Calculate the minimum cost among deletion, insertion, and substitution
            dp[i][j] = std::min({ 
                dp[i - 1][j] + 1,       // Deletion
                dp[i][j - 1] + 1,       // Insertion
                dp[i - 1][j - 1] + cost // Substitution
            });
        }
    }

    // Return the computed Levenshtein distance
    return dp[len1][len2];
}

// Function to process a chunk of records for a given long read and produce a graph object
Graph init_graph(const string& lr_name, std::vector<Record>& chunk,const string& lr_seq,bool verbose = false) {

    const int N = chunk.size();
    int lr_len = lr_seq.length();

    Graph G;
    std::map<string, boost::graph_traits<Graph>::vertex_descriptor> vertex_map;


    for (auto& node : chunk) {
        auto v = boost::add_vertex(boost::property<boost::vertex_index_t, string>(node.id), G);
        vertex_map[node.id] = v;
    }


    for(auto& i: chunk){
        
        for(auto& j: chunk){
            if( (i.id != j.id) && (i.start <= j.start && i.end < j.end) && (j.start <= i.end - MINOVERLAP + 1 ) && (j.end <= lr_len) ){
                std::vector<int> I = {j.start, i.end};
                std::vector<int> A = {I[0] - i.start, I[1] - I[0]}; // overlapping interval relative to i
                std::vector<int> B = {I[0] - j.start, I[1] - I[0]}; // overlapping interval relative to j
                

                std::string sub_i = i.sequence.substr(A[0], A[1]);
                std::string sub_j = j.sequence.substr(B[0], B[1]); 
                
                if (sub_i == sub_j) {

                    string j_overhang = j.sequence.substr(B[1], j.sequence.length() - B[1]); //portion of j's sequence which is not part of the overlap
                    string lr_subseq = lr_seq.substr(i.end, j_overhang.length());
                    int w = levenshteinDistance(j_overhang, lr_subseq);

                    boost::add_edge(vertex_map[i.id], vertex_map[j.id], boost::property<boost::edge_weight_t, int>(w), G);

                    if(verbose){

                        cout << "++++++++++++++++++++++++++++++++++\n" <<endl;
                        print_record(i,true);
                        print_record(j,true);
                        print_interval(I);
                        print_interval(A);
                        print_interval(B);
                        cout <<endl;

                        cout <<"j overhang ="<< j_overhang <<endl;
                        cout <<"lr substring ="<< lr_subseq <<endl;
                        cout <<"leven = " <<w<<endl;
                        cout << "\n\n";

                        cout << "LR:\n"<<lr_seq;
                    }

                }

                
            }
        }
    }

    return G;
    
}

void parseAndProcessData(const std::string& filepath,const std::unordered_map<std::string, std::string> LR_MAP) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filepath << std::endl;
        return;
    }
    

    std::string line;
    std::string currentPac = "";
    std::vector<Record> currentChunk;
    int mark = 5000;
    int count = 0;
    while (std::getline(file, line)) {
        
        std::istringstream iss(line);
        Record record;
        if (!(iss >> record.id >> record.pac >> record.start >> record.end >> record.sequence)) {
            std::cerr << "Error: Malformed line -> " << line << std::endl;
            continue;
        }
        if (record.pac == "*") {
            continue;
        }
        if (record.pac != currentPac) {
            if (!currentChunk.empty()) {
                init_graph(currentPac, currentChunk, LR_MAP.at(currentPac));
            }
            currentPac = record.pac;
            currentChunk.clear();
        }
        currentChunk.push_back(record);

        count ++;
        if(count >= mark){
            cout << "parsed "<< mark << " lines" <<endl;
            mark +=5000;
        }
    }

    if (!currentChunk.empty()) {
        Graph G = init_graph(currentPac, currentChunk, LR_MAP.at(currentPac));
    }

    file.close();
}


std::unordered_map<std::string, std::string> parseFasta(const std::string& filepath) {
    std::unordered_map<std::string, std::string> fastaDict;
    std::ifstream fastaFile(filepath);
    if (!fastaFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filepath << std::endl;
        return fastaDict;
    }
    

    std::string line, readName, sequence;
    while (std::getline(fastaFile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            
            if (!readName.empty()) {
                fastaDict[readName] = sequence;
            }
            readName = line.substr(1); // Remove '>'
            sequence.clear();         // Reset sequence for the new read
        } else {
            sequence += line; // Append the sequence line
        }
        
    }

    // Add the last read and sequence
    if (!readName.empty()) {
        fastaDict[readName] = sequence;
    }

    fastaFile.close();
    return fastaDict;
}



int main(int argc, char* argv[]) {

    
    if( argc != 3 ){
        cout << "Usage: ./build_graph <long_reads>.fasta <sl_raw_align.txt>"<< endl;
        exit(1);
    }
    cout << "Building Graphs...." << endl;
    std::unordered_map<std::string, std::string> LR_MAP = parseFasta(argv[1]);
    parseAndProcessData(argv[2] , LR_MAP);


    return 0;
}
