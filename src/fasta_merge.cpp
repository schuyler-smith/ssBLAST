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

std::vector<std::string> extract_keys(std::unordered_map<std::string, std::string> const& input_map) 
{
	std::vector<std::string> values;
	for (auto const& element : input_map) 
	{
  		values.push_back(element.first);
  	}
  return values;
}

std::vector<std::string> extract_values(std::unordered_map<std::string, std::string> const& input_map) 
{
	std::vector<std::string> values;
	for (auto const& element : input_map) 
	{
		values.push_back(element.second);
	}
	return values;
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::StringVector fasta_merge(
	Rcpp::StringVector fasta_file_path,
	const char * output_file
){
	Rcpp::Function directory("dir");
	Rcpp::Function file_path("normalizePath");
	Rcpp::Function basename("basename");
	Rcpp::Function paste("paste");

	Rcpp::StringVector 	read_IDs,
						read_seqs;
	std::string		fasta_line,
	 				read_ID,
					seq,
					line;

	std::unordered_map<std::string, std::string> reads;

	Rcpp::StringVector fasta_files;
	fasta_files = file_path(directory(fasta_file_path, Rcpp::Named("full.names", true)));
	if(fasta_files.length() == 0)
	{
		fasta_files = file_path(fasta_file_path);
	}
	
	size_t n_files = fasta_files.length();
	for(size_t i=0; i < n_files; ++i)
	{
		std::ifstream 	fasta_file(fasta_files[i]);
	    while(std::getline(fasta_file, line))
	    {
	    	std::stringstream   line_(line);
	    		std::getline(line_, fasta_line, '\t');
	    		if(fasta_line[0] == '>')
	    		{	  
	    			if(read_ID.compare("") != 0){ reads[read_ID] = seq; }
	    			read_ID = fasta_line;
	    			seq = "";
	    		} else {
	    			seq = seq + fasta_line;	    			
	    		}
	    }
	    reads[read_ID] = seq;
	}
	read_IDs = extract_keys(reads);	
	read_seqs = extract_values(reads);

	int lines = 0;
	remove(output_file);
	std::fstream output_file_path(output_file, std::ios::out | std::ios_base::app);
	std::stringstream output_line(std::ios_base::in | std::ios_base::out);
	size_t n_reads = read_IDs.length();
	for(size_t i=0; i < n_reads; ++i){
		output_line << read_IDs[i] << "\n" << read_seqs[i] << "\n" ;
		++lines;
		if(lines == 1000)
		{
			output_file_path << output_line.rdbuf();
			output_line.clear(); 
			lines = 0;
		}
	}
	output_file_path << output_line.rdbuf();
	output_file_path.close();

return(output_file);
}