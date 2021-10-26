/*
 *
 *	Author: Schuyler D. Smith
 *
 */

#include "BLAST_fastq_match.h"

#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <zlib.h>

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::CharacterVector aligned_BLAST_sequences(
    std::string FASTQ_file_path,
    Rcpp::StringVector aligned_reads,
    const char * output_path,
    bool aligned = false
){
    std::unordered_map<std::string, int> aligned_map;
    std::size_t num_aligned = aligned_reads.length();
    for(int i=0; i < num_aligned; ++i){
        aligned_map[Rcpp::as<std::string>(aligned_reads(i))] = 0;
    }
    
    bool compressed = {false};
    std::string ext = FASTQ_file_path.substr(FASTQ_file_path.rfind('.'));
    if(strcmp(ext.c_str(), ".gz") == 0){ compressed = true; } 
    else if(strcmp(ext.c_str(), ".bz2") == 0) { // bzip2 FASTQ_file_path
        Rcpp::Rcout << "bz file-type encription not supported.";
    } 
    if((strcmp(ext.c_str(), ".fastq") == 0) | (strcmp(ext.c_str(), ".fq") == 0) | compressed) {
        search_BLAST(FASTQ_file_path, aligned_map, output_path, aligned, compressed);
    }
    return(0);
}