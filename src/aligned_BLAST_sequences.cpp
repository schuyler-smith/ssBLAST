/*
 *
 *	Author: Schuyler D. Smith
 *
 */

#include <Rcpp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <zlib.h>

void search_aligned_BLAST(std::istream& FASTQ_file, std::unordered_map<std::string, int> aligned, const char * output_path)
{
	remove(output_path);
	std::string line, 
				seq_id;
	int lines_processed = 0, seq_line = 0;
	bool match = false;
	
	std::fstream output_file(output_path, std::ios::out | std::ios_base::app);
	std::stringstream output_line(std::ios_base::in | std::ios_base::out);
	while(std::getline(FASTQ_file, line))
	{
		if(seq_line == 4){ seq_line = 0; match = false; }
		if(seq_line > 0 && !(match)){ ++seq_line; continue; }
		if(seq_line == 0)
		{
			seq_id = line.substr(0, line.find(" ", 0)).erase(0,1);
			if(aligned.count(seq_id))
				{
					match = true;
				} else { ++seq_line; continue; }
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

void search_aligned_BLAST_gz(gzFile& FASTQ_file, std::unordered_map<std::string, int> aligned, const char * output_path)
{
	remove(output_path);
	std::string line, seq_id;
	int seq_line = 0;
	bool match = false;
 
	std::fstream output_file(output_path, std::ios::out | std::ios_base::app);
	static const unsigned BUFLEN = 1024;
	char buffer[BUFLEN];
	char* offset = buffer;
	for(;;)
	{
		int len = sizeof(buffer)-(offset-buffer);
		len = gzread(FASTQ_file, offset, len);
		if(len == 0) break;    
		char* cur = buffer;
		char* end = offset+len;
		for(char* eol; (cur < end) && (eol = std::find(cur, end, '\n')) < end; cur = eol + 1)
		{
			line = std::string(cur, eol);
			if(seq_line == 4){ seq_line = 0; match = false; }
			if(seq_line > 0 && !(match)){ ++seq_line; continue; }
			if(seq_line == 0)
			{
				seq_id = line.substr(0, line.find(" ", 0)).erase(0,1);
				if(aligned.count(seq_id))
  				{
  					match = true;
  				} else { ++seq_line; continue; }
			}
			output_file << line << "\n"; 
			++seq_line;
		}
		offset = std::copy(cur, end, buffer);
	}
	output_file << std::string(buffer, offset);
	output_file.close();
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
Rcpp::CharacterVector aligned_BLAST_sequences(
	std::string FASTQ_file_path,
	Rcpp::StringVector aligned,
	const char * output_path
){
	std::unordered_map<std::string, int> aligned_map;
	int num_aligned = aligned.length();
	for(int i=0; i < num_aligned; ++i)
    {
    	aligned_map[Rcpp::as<std::string>(aligned(i))] = 0;
    }
	std::string ext = FASTQ_file_path.substr(FASTQ_file_path.rfind('.'));
	if(strcmp(ext.c_str(), ".gz") == 0)
	{ 
		gzFile FASTQ_file = gzopen(FASTQ_file_path.c_str(), "rb");
		search_aligned_BLAST_gz(FASTQ_file, aligned_map, output_path);
		gzclose(FASTQ_file);
	} else if(strcmp(ext.c_str(), ".bz2") == 0) { // bzip2 FASTQ_file_path
		Rcpp::Rcout << "bz file-type encription not supported.";
	} else if((strcmp(ext.c_str(), ".fastq") == 0) | (strcmp(ext.c_str(), ".fq") == 0)) {
		std::ifstream FASTQ_file(FASTQ_file_path, std::ios::in);
		search_aligned_BLAST(FASTQ_file, aligned_map, output_path);
	}
	return(0);
}
