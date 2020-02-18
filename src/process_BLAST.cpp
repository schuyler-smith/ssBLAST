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

bool matches_cutoff(int x, int y){return(x>=y);}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::List process_BLAST(
    std::string BLAST_file_path,
    int min_length = 100,
    Rcpp::IntegerVector min_perc_ids = Rcpp::IntegerVector::create(100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90)
){
	std::string 	query,
					subject,
					e_value,
					line,
					line_,
					prior_query;

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

	Rcpp::StringVector BLAST_files;
	BLAST_files = file_path(directory(BLAST_file_path, Rcpp::Named("full.names", true)));
	if(BLAST_files.length() == 0){
		BLAST_files = file_path(BLAST_file_path);
	}
	
	size_t n_files = BLAST_files.length();
	std::vector<int> unique_reads(n_files);

	size_t n_cutoffs = min_perc_ids.length();
	Rcpp::List BLAST_counts(2 + n_cutoffs);
	BLAST_counts.fill(Rcpp::IntegerVector(n_files, 0));

	int min_perc_id = Rcpp::min(min_perc_ids);
	for(size_t i=0; i < n_files; ++i)
	{
		int unique_read = 0;
		std::vector<double> read_perc_id;

		std::ifstream 	BLAST_file(BLAST_files[i]);
	    while(std::getline(BLAST_file, line)) 
	    {
	    	std::stringstream   line_(line);
	    	std::getline(line_, query, '\t');
	    	line_ >> subject >> perc_id >> length >> mismatch >> gap >> 
	    		qstart >> qend >> sstart >> send >> e_value >> bitscore;
	    	if(length >= min_length)
		    {
			    if((query != prior_query) & (perc_id >= min_perc_id))
		   		{
		    		++unique_read;
		    		read_perc_id.push_back(perc_id);
		    	}
	    	}
	    }

	    unique_reads[i] = unique_read; BLAST_counts[1] = unique_reads;
	    for(size_t n=0; n < n_cutoffs; ++n){
	    	Rcpp::IntegerVector counts = BLAST_counts[2 + n];
	    	counts[i] = std::count_if(read_perc_id.begin(), 
							read_perc_id.end(), 
							[&](int const &x) {
									return(matches_cutoff(x, min_perc_ids[n]));
							});
	    	BLAST_counts[2+n] = counts;
	    }
	}
	BLAST_counts[0] = basename(BLAST_files);
	Rcpp::StringVector BLAST_counts_names(2 + n_cutoffs);
	for(size_t x = 0; x < n_cutoffs; ++x){BLAST_counts_names[2 + x] = std::to_string(min_perc_ids[x]);}
		BLAST_counts_names[0] = "Files"; BLAST_counts_names[1] = "Unique_Reads";
	BLAST_counts.attr("names") = BLAST_counts_names;
return(BLAST_counts);
}