#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>
using namespace boost;

using std::cout;
using std::endl;
using std::string;

const string EXAMPLE_LR = "pac_1";
const int MINOVERLAP = 10;
double bad_components = 0; //number of components where the shortest path was computed with a start position after the end position
double twice_aligned_short_reads = 0; //number of times a short read was aligned in more than one position to a single long read
double total_components = 0; //number of times a short read was aligned in more than one position to a single long read

using Graph = adjacency_list<
        vecS,                                    // Edge container (std::vector for adjacency list)
        vecS,                                    // Vertex container (std::vector for vertex list)
        undirectedS,                               // Directed graph
        property<vertex_index_t, string>, // Vertex properties: string id
        property<edge_weight_t, int>   // Edge properties integer weights (levenstien distance)
    >;

typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;

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


void write_graph_to_dot(const Graph& G, const std::map<std::string, graph_traits<Graph>::vertex_descriptor>& vertex_map, const std::string& filename) {
    std::ofstream dot_file(filename);

    dot_file << "digraph G {\n";

    // Set to track vertices that are part of edges
    std::set<graph_traits<Graph>::vertex_descriptor> connected_vertices;

    // Write edges and track connected vertices
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
        auto u = source(*ei, G); // Source vertex
        auto v = target(*ei, G); // Target vertex
        int weight = get(edge_weight, G, *ei); // Edge weight

        // Add vertices to the connected set
        connected_vertices.insert(u);
        connected_vertices.insert(v);

        // Write the edge
        dot_file << "\t" << u << " -> " << v << " [label=\"" << weight << "\"];" << std::endl;
    }

    // Write vertices that are part of edges
    for (const auto& pair : vertex_map) {
        if (connected_vertices.count(pair.second) > 0) { // Check if vertex is in the connected set
            dot_file << "\t" << pair.second << " [label=\"" << pair.first << " | " <<pair.second << "\"];" << std::endl;
        }
    }

    dot_file << "}\n";
    dot_file.close();

}



// Function to process a chunk of records for a given long read and produce a graph object, and a dictionary to map the vertex name to an id
std::tuple<Graph, std::map<string, graph_traits<Graph>::vertex_descriptor> > init_graph(

    const string& lr_name,
    std::vector<Record>& chunk,
    const string& lr_seq,
    bool verbose = false
    
    ) {

    const int N = chunk.size();
    int lr_len = lr_seq.length();

    Graph G;
    std::map<string, graph_traits<Graph>::vertex_descriptor> vertex_map;


    for (auto& node : chunk) {
        auto v = add_vertex(property<vertex_index_t, string>(node.id), G);
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

                    add_edge(vertex_map[i.id], vertex_map[j.id], property<edge_weight_t, int>(w), G);

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

    if(lr_name == EXAMPLE_LR){
        write_graph_to_dot(G ,vertex_map,"graph.dot");
    }

    return std::make_tuple(G,vertex_map);
    
}




string correct_read(
    Graph G,
    std::map<string, graph_traits<Graph>::vertex_descriptor> vertex_map,
    std::vector<Record>& chunk,
    const string& lr_name,
    const string& lr_seq,
    bool verbose = false
    
    ) {

    //map the vertex id (int) to the vertex name (string) found in the data 
    std::map<graph_traits<Graph>::vertex_descriptor, std::string> r_vertex_map;
    for (const auto& pair : vertex_map) {
        r_vertex_map[pair.second] = pair.first;  
    }
    // Create a vector to store the component index for each vertex
    std::vector<int> component(num_vertices(G));
    int num_components = connected_components(G, make_iterator_property_map(component.begin(), get(vertex_index, G)));
    // Map to store vertices for each component
    std::map<int, std::vector<Graph::vertex_descriptor>> component_vertices;
    for (auto v : make_iterator_range(vertices(G))) {
        component_vertices[component[v]].push_back(v);
    }
    // For each component, calculate the shortest path from the first node to the last node in the component list
    for (const auto& [component_id, vertices_list] : component_vertices) {
        total_components++;
        if (vertices_list.size() == 1) {
            continue;  // Skip components with only one vertex
        }
        // Extract the first and last vertices in the component list
        Vertex source = vertices_list.front();
        Vertex destination = vertices_list.back();

        // Define vectors for storing distances and predecessors
        std::vector<int> distance(num_vertices(G), std::numeric_limits<int>::max());
        std::vector<Vertex> predecessor(num_vertices(G), -1);

        // Run Dijkstra's algorithm
        dijkstra_shortest_paths(
            G,
            source,
            distance_map(
                make_iterator_property_map(distance.begin(), get(vertex_index, G))).predecessor_map(make_iterator_property_map(predecessor.begin(), get(vertex_index, G))
                )
            );


        std::vector<string> path;
        //get the id (first column in the alignment file) of each read in path
        for (graph_traits<Graph>::vertex_descriptor v = destination; v != source; v = predecessor[v]) {
            path.push_back(r_vertex_map[v]);  // Push the string ID
        }
        path.push_back(r_vertex_map[source]);
        std::vector<Record> subset;

        //get the node data from the path
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            for (const auto& record : chunk) {
                if (record.id == *it) {
                    subset.push_back(record);
                }
            }
        }
        
        //stop if a short read in the path 
        if(subset.size() != path.size()){
                
                twice_aligned_short_reads++;
                // one instance where a short read mapped to the same lr twice, meaning the reconstructed path contained 3 nodes instead of the expected 2.
                // id =ill.4339.1 | start = 692 | end =792
                // id =ill.8726.1 | start = 763 | end =863
                // id =ill.8726.1 | start = 984 | end =1084
                continue;
        }

        int l = subset.size();
        string s = subset[0].sequence;
        bool bad_component = false;
        for(int i = 0; i < l - 1; i++){
            int offset = subset[i].end - subset[i+1].start;
            if(offset >= 100){
                /*
                uncanny
                id =ill.2608.0 | start = 1006 | end =1106 | seq =  CTGACGTTCACGCTTACGTCCACACGGCATTCGGCAGATATTCCGCCGTATACGTTTGCCAGCGATGTGCAGGTTATGGTGATTAAGAAACAGGCGCTGG
                id =ill.3615.0 | start = 969 | end =1069 | seq =  TATTGTTGACATGCCAGCGGGTCGGGGAAACGTGATCCTGACGTTCACGCTTACGTCCACACGGCATTCGGCAGATATTCCGCCGTATACGTTTGCCAGC
                */
                bad_component = true;
                bad_components ++;
                s = "bad";
                break;
            }
            s += subset[i+1].sequence.substr(offset);
    
        }

        if(lr_name == EXAMPLE_LR and verbose){
                

                // Output the shortest distance and the path for the current component
                cout << "Component " << component_id << " - Shortest distance from " << path[path.size() - 1] << " to " << path[0] << " is: " << distance[destination] << endl;
                // Reconstruct the path from source to destination
                cout << "Path: \n";
                for(auto& it: subset){
                    cout << "\t";
                    print_record(it,true);
                }

                cout << "\ts = "<<s<<endl;

                cout << "_____________________________\n\n\n" << endl;
                

                // cout << "\ts = " << s;
                // cout << "\n------\n";



            }
        }

    return ""; 

    }

void correct_long_reads(const std::string& filepath,const std::unordered_map<std::string, std::string> LR_MAP) {
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
                string lr_seq = LR_MAP.at(currentPac);
                auto [G, vertex_map] = init_graph(currentPac, currentChunk, lr_seq);
                string s = correct_read(G,vertex_map, currentChunk, currentPac,lr_seq);
                

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
        string lr_seq = LR_MAP.at(currentPac);
        auto [G, vertex_map] = init_graph(currentPac, currentChunk, LR_MAP.at(currentPac));
        string s = correct_read(G,vertex_map, currentChunk, currentPac,lr_seq );

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
    correct_long_reads(argv[2] , LR_MAP);

    cout << "total components " <<total_components <<endl;
    cout << "bad components "<<bad_components<<" | " << 100* (bad_components/total_components)<< "%"<<endl;
    cout << "twice aligned short reads "<<twice_aligned_short_reads<<" | " << 100* (twice_aligned_short_reads/total_components)<< "%"<<endl;
    return 0;
}
