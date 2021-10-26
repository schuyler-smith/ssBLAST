#ifndef BLAST_FASTQ_MATCH_H_SDS
#define BLAST_FASTQ_MATCH_H_SDS

#include <Rcpp.h>
#include <string>
#include <zlib.h>

void search_BLAST(std::string, 
    std::unordered_map<std::string, int>, 
    const char *, 
    bool,
    bool
);

#endif /* BLAST_FASTQ_MATCH_H_SDS */