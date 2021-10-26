#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <zlib.h>

void search_BLAST(std::string FASTQ_file_path, 
    std::unordered_map<std::string, int> aligned_read_map, 
    const char * output_path, 
    bool aligned = false,
    bool compressed = false
){
    remove(output_path);
    std::string line;
    std::string seq_id;
    int lines_processed = {0}; 
    int seq_line = {0};
    bool match = {false};

    std::fstream output_file(output_path, std::ios::out | std::ios_base::app);
    if(!(compressed)){
        std::stringstream output_line(std::ios_base::in | std::ios_base::out);
        std::ifstream FASTQ_file(FASTQ_file_path, std::ios::in);
        while(std::getline(FASTQ_file, line)){
            if(seq_line == 4){ seq_line = 0; match = false; }
            if(seq_line == 0){
                seq_id = line.substr(0, line.find(" ", 0)).erase(0,1);
                if(aligned_read_map.count(seq_id)){ match = true; } 
            }
            ++seq_line;
            if(!(aligned) && match){ continue; }
            if(aligned && !(match)){ continue; }
            output_line << line << std::endl;
            ++lines_processed; 
            if(lines_processed == 1000){
                output_file << output_line.rdbuf();
                output_line.clear();
                lines_processed = 0;
            }
        }
        output_file << output_line.rdbuf();
    }
    if(compressed){
        static const unsigned BUFLEN = {1024};
        char buffer[BUFLEN];
        char* offset = buffer;
        gzFile FASTQ_file = gzopen(FASTQ_file_path.c_str(), "rb");
        for(;;){
            int len = sizeof(buffer)-(offset-buffer);
            len = gzread(FASTQ_file, offset, len);
            if(len == 0){ break; }    
            char* cur = buffer;
            char* end = offset + len;
            for(char* eol; (cur < end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1){
                line = std::string(cur, eol);
                if(seq_line == 4){ seq_line = 0; match = false; }
                if(seq_line == 0){
                    seq_id = line.substr(0, line.find(" ", 0)).erase(0,1);
                    if(aligned_read_map.count(seq_id)){ match = true; }
                }
                ++seq_line;
                if(!(aligned) && match){ continue; }
                if(aligned && !(match)){ continue; }
                output_file << line << std::endl; 
            }
            offset = std::copy(cur, end, buffer);
        }
        gzclose(FASTQ_file);
        output_file << std::string(buffer, offset);
    }
    output_file.close();
}