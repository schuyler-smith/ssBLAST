#' Extract reads not aligned in BLAST
#'
#' Compares fastq file with corresponding BLAST file to find reads that were not found in BLAST search and writes to new file.

#' @useDynLib ssBLAST
#' @usage unaligned_BLAST(fastq_file_path, BLAST_file_path, output_path = NULL)
#' @param fastq_file_path String that either represents a path to a single fastq file, or a directory containing multiple fastq files.
#' @param BLAST_file_path String that either represents a path to a single blast file, or a directory containing multiple blast files that correspond to `the fastq_file_path`.
#' @param output_path String for path that new files should be created. Default uses the `fastq_file_path`.
#' @return fastq
#' @import Rcpp
#' @importFrom data.table fread
#' @export
#' 

# Sys.setenv("PKG_LIBS"="-lboost_iostreams")
unaligned_BLAST <- function(fastq_file_path, BLAST_file_path, output_path = NULL){
  
  if(length(fastq_file_path) == 1 && dir.exists(fastq_file_path)){ 
    fastq_file_path <- gsub('/$','',fastq_file_path)
    fastq_file <- vector()
    for(extension in c(".fastq.gz$", ".fastq.bz2$", ".fastq$", ".fq$")){
      fastq_file <- append(fastq_file, dir(fastq_file_path, extension, full.names = TRUE))
    }
  } else {fastq_file <- fastq_file_path}
  fastq_file <- sort(normalizePath(fastq_file))
  if(length(BLAST_file_path) == 1 && dir.exists(BLAST_file_path)){ 
    BLAST_file_path <- gsub('/$','',BLAST_file_path)
    BLAST_file <- dir(BLAST_file_path, full.names = TRUE)
  } else {BLAST_file <- BLAST_file_path}
  BLAST_file <- sort(normalizePath(BLAST_file))
  if(length(fastq_file) > 1 || length(BLAST_file) > 1){
    if(length(fastq_file) != length(BLAST_file)){
      warning('Number of submitted FASTQ and BLAST files are not equal\nonly files with corresponding names will be used.')
      fastqs <- basename(sapply(strsplit(fastq_file,'\\.'),"[[", 1))
      blasts <- basename(sapply(strsplit(BLAST_file,'\\.'),"[[", 1))
      fastq_names <- unlist(sapply(fastqs, FUN = function(x){ grep(x, blasts) }))
      blast_names <- unlist(sapply(blasts, FUN = function(x){ grep(x, fastqs) }))
      if(length(fastq_names) > length(blast_names)){
        fastq_file <- fastq_file[which(fastqs %in% names(fastq_names))]
        BLAST_file <- BLAST_file[fastq_names]
      } else {
        BLAST_file <- BLAST_file[which(blasts %in% names(blast_names))]
        fastq_file <- fastq_file[blast_names]
      }
      if(length(fastq_file) != length(BLAST_file) || length(BLAST_file) == 0 || length(fastq_file) == 0){
        stop('All fastq files should have a corresponding BLAST file.')
      }
    }
  }
  if(!is.null(output_path)){
    output_path <- normalizePath(output_path)
    if(dir.exists(output_path)){
      for(file in seq_along(fastq_file)){
        unaligned_BLAST_sequences(fastq_file[file], 
              unique(data.table::fread(BLAST_file[file])[[1]]), 
              file.path(gsub('/$','',normalizePath(output_path)), paste("unaligned_", 
                   basename(sapply(strsplit(fastq_file[file],'\\.'),"[[", 1)), 
                   ".fastq", sep = '')))
      } 
    } else {
      for(file in seq_along(fastq_file)){
        unaligned_BLAST_sequences(fastq_file[file], 
              unique(data.table::fread(BLAST_file[file])[[1]]), 
              output_path)
      } 
    }
  } else {
    for(file in seq_along(fastq_file)){
      unaligned_BLAST_sequences(fastq_file[file], 
        unique(data.table::fread(BLAST_file[file])[[1]]), 
        file.path(dirname(fastq_file[file]), paste("unaligned_", 
                                                   basename(sapply(strsplit(fastq_file[file],'\\.'),"[[", 1)), 
                                                   ".fastq", sep = '')))
    } 
  }
}

