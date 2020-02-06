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

std::vector<std::string> extract_keys(std::unordered_map<std::string, double> const& input_map) 
{
	std::vector<std::string> values;
	for (auto const& element : input_map) 
	{
  		values.push_back(element.first);
  	}
  return values;
}

std::vector<double> extract_values(std::unordered_map<std::string, double> const& input_map) 
{
	std::vector<double> values;
	for (auto const& element : input_map) 
	{
		values.push_back(element.second);
	}
	return values;
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::DataFrame BLAST_db_match(
	std::string BLAST_file_path,
	std::string db_file_path,
	int min_length = 100,
	int min_id = 98
){

	std::string 	query,
					prior_query,
					db_match,
					e_value,
					line,
					line_,
					db_ID,
					not_db_ID;;

	int 		length,
					mismatch,
					gap,
					qstart,
					qend,
					sstart,
					send,
					count;
	double 	perc_id,
					bitscore,
					prior_perc_id,
					prior_bitscore;
	std::unordered_map<std::string, double> dbs;
	
	Rcpp::Function directory("dir");
	Rcpp::Function file_path("normalizePath");
	Rcpp::Function basename("basename");
	
	Rcpp::StringVector db_files;
	db_files = file_path(db_file_path);
	std::ifstream 	db_file(db_files[0]);
    while(std::getline(db_file, line))
    {
    	std::stringstream   line_(line);
    	std::getline(line_, db_ID, '\t');
    	line_ >> not_db_ID;
    	if(db_ID[0] == '>')
    	{
    		dbs[db_ID.erase(0, 1)] = 0;
    	}
    }

	Rcpp::StringVector BLAST_files = directory(BLAST_file_path);
	Rcpp::StringVector namevec;
	if(BLAST_files.length() > 0)
	{
		BLAST_files = file_path(directory(BLAST_file_path, Rcpp::Named("full.names", true)));
		namevec = directory(BLAST_file_path);
	} else {
		BLAST_files = file_path(BLAST_file_path);
		namevec = basename(BLAST_file_path);
	}
	namevec.push_front("db");
	size_t n_files = BLAST_files.length();

	Rcpp::List db_counts(1 + n_files);
	db_counts[0] = extract_keys(dbs);

	for(size_t i=0; i < n_files; ++i)
	{
		std::unordered_map<std::string, double> found_dbs = dbs;
		std::ifstream 	BLAST_file(BLAST_files[i]);
	    while(std::getline(BLAST_file, line))
	    {
	    	std::stringstream   line_(line);
	    	std::getline(line_, query, '\t');
	    	line_ >> db_match >> perc_id >> length >> mismatch >> gap >>
	    		qstart >> qend >> sstart >> send >> e_value >> bitscore;
	    	if((query != prior_query) &
	    		(length >= min_length) &
	    		(perc_id >= min_id))
	    	{
  				if(found_dbs.count(db_match))
  				{
  					found_dbs[db_match]++;
  				}
	    	}
	    	prior_query = query;
	    }
    	db_counts[i+1] = extract_values(found_dbs);
	}
	db_counts.attr("names") = namevec;

	Rcpp::DataFrame return_dataframe(db_counts);
return(return_dataframe);
}

