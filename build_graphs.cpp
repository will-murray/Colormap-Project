#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
using std::cout;
using std::endl;
using std::string;

const int MINOVERLAP = 10;


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
    cout<< "["<<I[0] << ", "<<I[1] <<"]";
}
// Function to process a chunk of records for a given long read
void processChunk(const std::string& pac,  std::vector<Record>& chunk, const std::unordered_map<std::string, std::string> LR_MAP) {
    
    for(auto& i: chunk){
        
        for(auto& j: chunk){
            if( (i.id != j.id) && (i.start <= j.start && i.end < j.end) && (j.start <= i.end - MINOVERLAP + 1 ) ){
                std::vector<int> I = {j.start, i.end};
                std::vector<int> A = {I[0] - i.start, I[1] - I[0]}; // overlapping interval relative to i
                std::vector<int> B = {I[0] - j.start, I[1] - I[0]}; // overlapping interval relative to j
                

                std::string sub_i = i.sequence.substr(A[0], A[1]);
                std::string sub_j = j.sequence.substr(B[0], B[1]); 
                if (sub_i == sub_j) {
                    
                }

                

            }
        }
    }
}

void parseAndProcessData(const std::string& filepath,const std::unordered_map<std::string, std::string> LR_MAP) {
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filepath << std::endl;
        return;
    }
    cout << "opened " <<filepath << " succefully..."<<endl;

    std::string line;
    std::string currentPac = "";
    std::vector<Record> currentChunk;

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
                processChunk(currentPac, currentChunk, LR_MAP);
            }
            currentPac = record.pac;
            currentChunk.clear();
        }
        currentChunk.push_back(record);
    }

    if (!currentChunk.empty()) {
        processChunk(currentPac, currentChunk, LR_MAP);
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
    std::unordered_map<std::string, std::string> LR_MAP = parseFasta(argv[1]);
    parseAndProcessData(argv[2] , LR_MAP);

    return 0;
}
