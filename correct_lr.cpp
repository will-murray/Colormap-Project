#include <iostream>
#include <fstream>
#include <regex>

void correct(std::string data_fname,std::string freq_fname){
    std::ifstream data_file(data_fname); // Open the alignment file
    std::ifstream freq_file(data_fname); // Open the file

    std::string line;
    std::regex regex(R"(\s(\S+))"); //matches the lr_name
    std::string lr_name;
    lr_name = "*";
    int tokens;

    if( data_file.is_open() ){
        while( std::getline( data_file, line ) )
        {
            std::smatch match;
            std::regex_search( line, match, regex );
            if( match[1] != lr_name ){
                lr_name = match[1];
            } 

        }
        

    }

    


    return;
}


//compliation: g++ correct_lr.cpp -o correct_lr
int main(int argc, char* argv[]) {
    std::string freq_file;
    std::string data_file;
    std::string lr_file;

    if (argc == 1) {
        // Test mode
        freq_file = "test_data/freq.txt";
        data_file = "test_data/sl_raw_align.txt";
        lr_file = "test_data/pac.fasta";
    } else if (argc == 4) {
        freq_file = argv[1];
        data_file = argv[2];
        lr_file = argv[3];   
    } else{
        std::cout<<"wrong usage BOZO"<<std::endl;exit(1);
    }


    correct(data_file,freq_file);
    return 0;


}

