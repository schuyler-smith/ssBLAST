// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BLAST_db_match
Rcpp::DataFrame BLAST_db_match(std::string BLAST_file_path, std::string db_file_path, int min_length, int min_id);
RcppExport SEXP _ssBLAST_BLAST_db_match(SEXP BLAST_file_pathSEXP, SEXP db_file_pathSEXP, SEXP min_lengthSEXP, SEXP min_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type BLAST_file_path(BLAST_file_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type db_file_path(db_file_pathSEXP);
    Rcpp::traits::input_parameter< int >::type min_length(min_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type min_id(min_idSEXP);
    rcpp_result_gen = Rcpp::wrap(BLAST_db_match(BLAST_file_path, db_file_path, min_length, min_id));
    return rcpp_result_gen;
END_RCPP
}
// aligned_BLAST_sequences
Rcpp::CharacterVector aligned_BLAST_sequences(std::string FASTQ_file_path, Rcpp::StringVector aligned_reads, const char * output_path, bool aligned);
RcppExport SEXP _ssBLAST_aligned_BLAST_sequences(SEXP FASTQ_file_pathSEXP, SEXP aligned_readsSEXP, SEXP output_pathSEXP, SEXP alignedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type FASTQ_file_path(FASTQ_file_pathSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type aligned_reads(aligned_readsSEXP);
    Rcpp::traits::input_parameter< const char * >::type output_path(output_pathSEXP);
    Rcpp::traits::input_parameter< bool >::type aligned(alignedSEXP);
    rcpp_result_gen = Rcpp::wrap(aligned_BLAST_sequences(FASTQ_file_path, aligned_reads, output_path, aligned));
    return rcpp_result_gen;
END_RCPP
}
// derep_sequences
int derep_sequences(std::string FASTQ_file_path, std::string output_path, bool fasta, int n);
RcppExport SEXP _ssBLAST_derep_sequences(SEXP FASTQ_file_pathSEXP, SEXP output_pathSEXP, SEXP fastaSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type FASTQ_file_path(FASTQ_file_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_path(output_pathSEXP);
    Rcpp::traits::input_parameter< bool >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(derep_sequences(FASTQ_file_path, output_path, fasta, n));
    return rcpp_result_gen;
END_RCPP
}
// fasta_merge
Rcpp::StringVector fasta_merge(Rcpp::StringVector fasta_file_path, const char * output_file);
RcppExport SEXP _ssBLAST_fasta_merge(SEXP fasta_file_pathSEXP, SEXP output_fileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type fasta_file_path(fasta_file_pathSEXP);
    Rcpp::traits::input_parameter< const char * >::type output_file(output_fileSEXP);
    rcpp_result_gen = Rcpp::wrap(fasta_merge(fasta_file_path, output_file));
    return rcpp_result_gen;
END_RCPP
}
// fasta_seq_names
Rcpp::List fasta_seq_names(Rcpp::StringVector fasta_file_path);
RcppExport SEXP _ssBLAST_fasta_seq_names(SEXP fasta_file_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type fasta_file_path(fasta_file_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(fasta_seq_names(fasta_file_path));
    return rcpp_result_gen;
END_RCPP
}
// fasta_seq_reads
Rcpp::List fasta_seq_reads(Rcpp::StringVector fasta_file_path);
RcppExport SEXP _ssBLAST_fasta_seq_reads(SEXP fasta_file_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type fasta_file_path(fasta_file_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(fasta_seq_reads(fasta_file_path));
    return rcpp_result_gen;
END_RCPP
}
// filter_blast
Rcpp::StringVector filter_blast(const std::string& BLAST_file_path, const std::string& target_query, const std::string& target_subject, double perc_id_cutoff, int min_length, int max_qstart, int min_qend, int max_sstart, int min_send);
RcppExport SEXP _ssBLAST_filter_blast(SEXP BLAST_file_pathSEXP, SEXP target_querySEXP, SEXP target_subjectSEXP, SEXP perc_id_cutoffSEXP, SEXP min_lengthSEXP, SEXP max_qstartSEXP, SEXP min_qendSEXP, SEXP max_sstartSEXP, SEXP min_sendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type BLAST_file_path(BLAST_file_pathSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type target_query(target_querySEXP);
    Rcpp::traits::input_parameter< const std::string& >::type target_subject(target_subjectSEXP);
    Rcpp::traits::input_parameter< double >::type perc_id_cutoff(perc_id_cutoffSEXP);
    Rcpp::traits::input_parameter< int >::type min_length(min_lengthSEXP);
    Rcpp::traits::input_parameter< int >::type max_qstart(max_qstartSEXP);
    Rcpp::traits::input_parameter< int >::type min_qend(min_qendSEXP);
    Rcpp::traits::input_parameter< int >::type max_sstart(max_sstartSEXP);
    Rcpp::traits::input_parameter< int >::type min_send(min_sendSEXP);
    rcpp_result_gen = Rcpp::wrap(filter_blast(BLAST_file_path, target_query, target_subject, perc_id_cutoff, min_length, max_qstart, min_qend, max_sstart, min_send));
    return rcpp_result_gen;
END_RCPP
}
// process_BLAST
Rcpp::List process_BLAST(std::string BLAST_file_path, int min_length, Rcpp::IntegerVector min_perc_ids);
RcppExport SEXP _ssBLAST_process_BLAST(SEXP BLAST_file_pathSEXP, SEXP min_lengthSEXP, SEXP min_perc_idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type BLAST_file_path(BLAST_file_pathSEXP);
    Rcpp::traits::input_parameter< int >::type min_length(min_lengthSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type min_perc_ids(min_perc_idsSEXP);
    rcpp_result_gen = Rcpp::wrap(process_BLAST(BLAST_file_path, min_length, min_perc_ids));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ssBLAST_BLAST_db_match", (DL_FUNC) &_ssBLAST_BLAST_db_match, 4},
    {"_ssBLAST_aligned_BLAST_sequences", (DL_FUNC) &_ssBLAST_aligned_BLAST_sequences, 4},
    {"_ssBLAST_derep_sequences", (DL_FUNC) &_ssBLAST_derep_sequences, 4},
    {"_ssBLAST_fasta_merge", (DL_FUNC) &_ssBLAST_fasta_merge, 2},
    {"_ssBLAST_fasta_seq_names", (DL_FUNC) &_ssBLAST_fasta_seq_names, 1},
    {"_ssBLAST_fasta_seq_reads", (DL_FUNC) &_ssBLAST_fasta_seq_reads, 1},
    {"_ssBLAST_filter_blast", (DL_FUNC) &_ssBLAST_filter_blast, 9},
    {"_ssBLAST_process_BLAST", (DL_FUNC) &_ssBLAST_process_BLAST, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_ssBLAST(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
