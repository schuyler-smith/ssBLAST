/*
 *
 *	Author: Schuyler D. Smith
 *
 */

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

//' @author Schuyler D. Smith
// [[Rcpp::export]]

Rcpp::StringVector parse_blast(
	const std::string& BLAST_file_path,
	const std::string& target_query = "",
	const std::string& target_subject = "",
	double perc_id_cutoff = 90,
	int min_length = 100,
	int max_qstart = -1,
	int min_qend = -1,
	int max_sstart = -1,
	int min_send = -1
){
	std::string 	subject,
					e_value;

	int 			length,
					mismatch,
					gap,
					qstart,
					qend,
					sstart,
					send;

	double 			perc_id,
					bitscore;

	Rcpp::Function directory("dir");
	Rcpp::Function file_path("normalizePath");
	Rcpp::Function basename("basename");

	std::ifstream 	file("/home/schuyler/DARTE/BLAST/resfinder/Q3-X1P66N3BFD7_S76_L001_R2_001.blast");
	std::string 	line;

	Rcpp::StringVector BLAST_files = directory(BLAST_file_path);
	Rcpp::StringVector namevec;
	if(BLAST_files.length() > 0)
	{
		BLAST_files = file_path(directory(BLAST_file_path, Rcpp::Named("full.names", true)));
	} else {
		BLAST_files = file_path(BLAST_file_path);
	}

	std::string output_path;
	size_t n_files = BLAST_files.length();
	for(size_t i=0; i < n_files; ++i)
	{
		output_path = BLAST_files[i];
		output_path = output_path + ".filtered";
		std::fstream output_file(output_path, std::ios::out | std::ios_base::app);
		std::stringstream output_line(std::ios_base::in | std::ios_base::out);
		std::ifstream 	BLAST_file(BLAST_files[i]);
	    while(std::getline(BLAST_file, line))
	    {
			std::stringstream   line_s(line);
			std::string 		query;
			std::getline(line_s, query, '\t');
			line_s >> subject >> perc_id >> length >> mismatch >> gap >> 
			qstart >> qend >> sstart >> send >> e_value >> bitscore;
			if(target_query.compare("") != 0){
				if(query.compare(target_query) != 0){ continue; }
			}
			if(target_subject.compare("") != 0){
				if(subject.compare(target_subject) != 0){ continue; }
			}
			if(!(perc_id >= perc_id_cutoff)){ continue; } 
			if(!(length >= min_length)){ continue; } 
			if(max_qstart != -1){
				if(qstart > max_qstart){ continue; }
			}
			if(min_qend != -1){
				if(qend < min_qend){ continue; }
			}
			if(max_sstart != -1){
				if(sstart > max_sstart){ continue; }
			}
			if(min_send != -1){
				Rcpp::Rcout << line;
				if(send < min_send){ continue; }
			}
			output_line << line << '\n' ;
			output_file << output_line.rdbuf();
			output_line.clear();
		}
		output_file.close();
	}
return(0);
}