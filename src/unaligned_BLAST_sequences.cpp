/*
 *
 *	Author: Schuyler D. Smith
 *
 */

#include <Rcpp.h>
#include <string>

#include <fstream>
#include <sstream>
#include "filtering_stream.hpp"
#include "gzip.hpp"
#include "bzip2.hpp"
// #include <boost/iostreams/copy.hpp>

void search_BLAST(std::istream& FASTQ_file, Rcpp::StringVector aligned, const char * output_path)
{
	remove(output_path);
	std::string line,
				seq_id;
	int lines_processed = 0,
		seq_line = 0, 
		BLAST_seq = 0;
	bool match = false;

	auto output_file = std::fstream(output_path, std::ios::out | std::ios_base::app);
	std::stringstream output_line(std::ios_base::in | std::ios_base::out);
	while(std::getline(FASTQ_file, line)) 
	{
		if(seq_line == 4){ seq_line = 0; match = false; }
		if(seq_line > 0 && match){ ++seq_line; continue; }
		if(seq_line == 0)
		{
			seq_id = line.substr(0, line.find(" ", 0));
			if(strcmp(seq_id.erase(0,1).c_str(), aligned[BLAST_seq]) == 0)
			{
				match = true; ++BLAST_seq; ++seq_line; continue;
			}
		}
		output_line << line << "\n";
		++lines_processed; ++seq_line;
		if(lines_processed == 1000)
		{
			output_file << output_line.rdbuf();
			output_line.clear(); lines_processed = 0;
		}
	}
	output_file << output_line.rdbuf();
	output_file.close();
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::CharacterVector unaligned_BLAST_sequences(
	std::string FASTQ_file_path,
	Rcpp::StringVector aligned,
	const char * output_path
){
	std::string ext = FASTQ_file_path.substr(FASTQ_file_path.rfind('.'));
	if(strcmp(ext.c_str(), ".gz")==0 || strcmp(ext.c_str(), ".bz2")==0) { 
		std::ifstream file(FASTQ_file_path.c_str(), std::ios_base::in | std::ios_base::binary);
		if (!file) {Rcpp::Rcout << "Error reading compressed fastq file: " << FASTQ_file_path << "\n"; return R_NilValue;}
		boost::iostreams::filtering_istream FASTQ_file;
		if (strcmp(ext.c_str(), ".gz")==0) 
		{ // gzip FASTQ_file_path
			FASTQ_file.push(boost::iostreams::gzip_decompressor());
		} else { // bzip2 FASTQ_file_path
			FASTQ_file.push(boost::iostreams::bzip2_decompressor());
		}
		FASTQ_file.push(file);
		if(!FASTQ_file) {Rcpp::Rcout << "Error reading fastq file: " << FASTQ_file_path << "\n"; return R_NilValue;}
		search_BLAST(FASTQ_file, aligned, output_path);
	} else {
		std::ifstream FASTQ_file(FASTQ_file_path);
		search_BLAST(FASTQ_file, aligned, output_path);
	}
	
	return(0);
}
