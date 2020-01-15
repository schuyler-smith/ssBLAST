#' Dereplicate reads in FASTQ file.
#'
#' Dereplicates read from a FASTQ file. Keeps quality score according to the first encounter of the sequence, and 
#' sequence ID from the corresponding read that owns the quality score.

#' @useDynLib ssBLAST
#' @usage dereplicate_sequences(fastq_file_path, output_path = NULL, fasta, n)
#' @param fastq_file_path String that either represents a path to a single fastq file, or a directory containing multiple fastq files.
#' @param output_path String for path that new files should be created. Default uses the `fastq_file_path`.
#' @param fasta if `TRUE` then the sequences are output in fasta format, fastq if `FALSE`.
#' @param n number of sequences to return (most abundant). 0 (default) returns all reads.
#' @return fastq
#' @import Rcpp
#' @export
#' 

dereplicate_sequences <- function(fastq_file_path, output_path = NULL, fasta = FALSE, n = 0){
  
  if(length(fastq_file_path) == 1 && dir.exists(fastq_file_path)){ 
    fastq_file_path <- gsub('/$','',fastq_file_path)
    fastq_file <- vector()
    for(extension in c(".fastq.gz$", ".fastq.bz2$", ".fastq$")){
      fastq_file <- append(fastq_file, dir(fastq_file_path, extension, full.names = TRUE))
    }
  } else {fastq_file <- fastq_file_path}
  fastq_file <- normalizePath(fastq_file)
  
  if(!is.null(output_path)){
    output_path <- normalizePath(output_path)
    if(dir.exists(output_path)){
      output_file <- file.path(gsub('/$','',normalizePath(output_path)), paste("dereplicated_", 
                                                                basename(sapply(strsplit(fastq_file,'\\.'),"[[", 1)), 
                                                                ".fastq", sep = ''))
    } 
    } else {
      output_file <- file.path(dirname(fastq_file), paste("dereplicated_", 
                                                                   basename(sapply(strsplit(fastq_file,'\\.'),"[[", 1)), 
                                                                   ".fastq", sep = ''))
    }
    for(file in seq_along(fastq_file)){
      derep_sequences(fastq_file[file], 
                      output_file[file],
                      fasta,
                      n)
    } 

}
  