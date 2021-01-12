/*
 *
 *	Author: Schuyler D. Smith
 *
 */

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::List fasta_seq_names(
	Rcpp::StringVector fasta_file_path
){
	std::string 	seq_ID,
					not_seq_ID,
					line;

	Rcpp::Function directory("dir");
	Rcpp::Function file_path("normalizePath");
	Rcpp::Function basename("basename");

	Rcpp::StringVector fasta_files;
	Rcpp::StringVector namevec;

	fasta_files = file_path(directory(fasta_file_path, Rcpp::Named("full.names", true)));
	if(fasta_files.length() > 0)
	{ 
		namevec = directory(fasta_file_path);
	} else {
		fasta_files = file_path(fasta_file_path);
		namevec = basename(fasta_file_path);
	}

	size_t n_files = fasta_files.length();
	Rcpp::List fasta_seqs(n_files);
	for(size_t i=0; i < n_files; ++i)
	{
		Rcpp::StringVector seq_names;
		std::ifstream fasta_file(fasta_files[i]);
	    while(std::getline(fasta_file, line))
	    {
	    	std::stringstream   line_(line);
	    	std::getline(line_, seq_ID, '\t');
	    	line_ >> not_seq_ID;
	    	if(seq_ID[0] == '>')
	    	{
	    		seq_names.push_back(seq_ID.erase(0, 1));
	    	}
	    }
    	fasta_seqs[i] = seq_names;
	}
	fasta_seqs.attr("names") = namevec;
return(fasta_seqs);
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::List fasta_seq_reads(
	Rcpp::StringVector fasta_file_path
){
	Rcpp::Function directory("dir");
	Rcpp::Function file_path("normalizePath");
	Rcpp::Function basename("basename");

	Rcpp::StringVector fasta_files;
	Rcpp::StringVector namevec;

	fasta_files = file_path(directory(fasta_file_path, Rcpp::Named("full.names", true)));
	if(fasta_files.length() > 0)
	{ 
		namevec = directory(fasta_file_path);
	} else {
		fasta_files = file_path(fasta_file_path);
		namevec = basename(fasta_file_path);
	}

	size_t n_files = fasta_files.length();
	Rcpp::List fasta_seqs(n_files);
	for(size_t i=0; i < n_files; ++i)
	{	
		std::string 	read_seq,
						read_line,
						line;
		bool file_start(true);
		Rcpp::StringVector read_seqs;
		std::ifstream fasta_file(fasta_files[i]);
	    while(std::getline(fasta_file, line))
	    {
	    	std::stringstream line_(line);
	    	std::getline(line_, read_line, '\t');
	    	if(read_line[0] == '>')
	    	{
	    		if(!file_start)
	    		{
	    			read_seqs.push_back(read_seq);
	    			read_seq.clear();
	    		} else {
		    		file_start = false;	
	    		}
	    	} else {
	    		read_seq.append(read_line);
	    	}
	    }
	    read_seqs.push_back(read_seq);
    	fasta_seqs[i] = read_seqs;
	}
	fasta_seqs.attr("names") = namevec;
return(fasta_seqs);
}