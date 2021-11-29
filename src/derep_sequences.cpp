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

bool desc_sort(std::pair<std::string, int> a, std::pair<std::string, int> b) {
    return a.second > b.second;
}

void process_seqs(std::istream& FASTQ_file, const char * output_path, bool fasta, int n)
{
	remove(output_path);
	std::string line, seq_id, seq;
	int seq_line = 1;
	bool match = false;
	std::unordered_map<std::string, int> seq_counts;
	std::unordered_map<std::string, std::string> seq_ids;
	std::unordered_map<std::string, std::string> q_scores;

	while(std::getline(FASTQ_file, line)) 
	{
		if(seq_line == 1){ seq_id = line;}
		if(seq_line == 2)
		{
			seq = line;
			if(seq_counts.count(seq))
			{
				seq_counts[seq]++;
			} else {
				seq_counts[seq] = 1;
				seq_ids[seq] = seq_id;
				match = true;
			}
		}
		if(seq_line == 4)
		{ 
			if(match){ q_scores[seq] = line; } 
			seq_line = 1; 
			match = false; 
			continue; 
		}
		++seq_line;
	}
	std::vector<std::pair<std::string, int>> sorted_seq_counts(seq_counts.begin(), seq_counts.end());
	std::sort(sorted_seq_counts.begin(), sorted_seq_counts.end(), desc_sort);

	std::fstream output_file(output_path, std::ios::out | std::ios_base::app);
	if(n == 0){ n = sorted_seq_counts.size(); }
	for (int i = 0; i < n; ++i)
	{
		if(fasta){
			output_file << ">" << seq_ids[sorted_seq_counts[i].first].erase(0,1) << "\t" << sorted_seq_counts[i].second << "\n";
			output_file << sorted_seq_counts[i].first << std::endl;
		} else {
			output_file << seq_ids[sorted_seq_counts[i].first] << "\t" << sorted_seq_counts[i].second << "\n";
			output_file << sorted_seq_counts[i].first << "\n" << "+" << "\n";
			output_file << q_scores[sorted_seq_counts[i].first] << "\n";			
		}
	}
	output_file.close();
}

void process_seqs_gz(gzFile& FASTQ_file, const char * output_path, bool fasta, int n)
{
	remove(output_path);
	std::string line, seq_id, seq;
	int seq_line = 1;
	bool match = false;
	std::unordered_map<std::string, int> seq_counts;
	std::unordered_map<std::string, std::string> seq_ids;
	std::unordered_map<std::string, std::string> q_scores;

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
			if(seq_line == 1){ seq_id = line;}
			if(seq_line == 2)
			{
				seq = line;
				if(seq_counts.count(seq))
				{
					seq_counts[seq]++;
				} else {
					seq_counts[seq] = 1;
					seq_ids[seq] = seq_id;
					match = true;
				}
			}
			if(seq_line == 4)
			{ 
				if(match){ q_scores[seq] = line; } 
				seq_line = 1; 
				match = false; 
				continue; 
			}
			++seq_line;
		}
		offset = std::copy(cur, end, buffer);
	}
	std::vector<std::pair<std::string, int>> sorted_seq_counts(seq_counts.begin(), seq_counts.end());
	std::sort(sorted_seq_counts.begin(), sorted_seq_counts.end(), desc_sort);

	std::fstream output_file(output_path, std::ios::out | std::ios_base::app);
	if(n == 0){ n = sorted_seq_counts.size(); }
	for (int i = 0; i < n; ++i)
	{
		if(fasta){
			output_file << ">" << seq_ids[sorted_seq_counts[i].first].erase(0,1) << "\t" << sorted_seq_counts[i].second << "\n";
			output_file << sorted_seq_counts[i].first << std::endl;
		} else {
			output_file << seq_ids[sorted_seq_counts[i].first] << "\t" << sorted_seq_counts[i].second << "\n";
			output_file << sorted_seq_counts[i].first << "\n" << "+" "\n";
			output_file << q_scores[sorted_seq_counts[i].first] << std::endl;			
		}
	}
	output_file.close();
}

//' @author Schuyler D. Smith
// [[Rcpp::export]]
int derep_sequences(
	std::string FASTQ_file_path,
	std::string output_path, 
	bool fasta,
	int n
){
	const char * output_file;
	if(fasta){ output_path = output_path.substr(0, output_path.find(".")) + ".fasta"; }
	output_file = output_path.c_str();

	std::string ext = FASTQ_file_path.substr(FASTQ_file_path.rfind('.'));
	if(strcmp(ext.c_str(), ".gz") == 0) { 
		gzFile FASTQ_file = gzopen(FASTQ_file_path.c_str(), "rb");
		process_seqs_gz(FASTQ_file, output_file, fasta, n);
		gzclose(FASTQ_file);
	} else if(strcmp(ext.c_str(), ".bz2") == 0) {
		Rcpp::Rcout << "bz file-type encription not supported.";
	} else if(strcmp(ext.c_str(), ".fastq") == 0) {
		std::ifstream FASTQ_file(FASTQ_file_path, std::ios::in);
		process_seqs(FASTQ_file, output_file, fasta, n);
	}
	return(0);
}