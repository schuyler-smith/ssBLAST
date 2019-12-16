#' Extract reads not aligned in BLAST
#'
#' Compares fastq file with corresponding BLAST file to find reads that were not found in BLAST search and writes to new file.

#' @useDynLib ssBLAST
#' @usage unaligned_BLAST(fastq_file_path, BLAST_file_path, output_path = NULL)
#' @param fastq_file_path String that either represents a path to a single fastq file, or a directory containing multiple fastq files.
#' @param BLAST_file_path String that either represents a path to a single blast file, or a directory containing multiple blast files that correspond to `the fastq_file_path`.
#' @param output_path String for path that new files should be created. Default uses the `fastq_file_path`.
#' @return fastq
#' @importFrom data.table fread
#' @export
#' 

# Sys.setenv("PKG_LIBS"="-lboost_iostreams")
unaligned_BLAST <- function(fastq_file_path, BLAST_file_path, output_path = NULL){
  
  if(length(fastq_file_path) == 1 && dir.exists(fastq_file_path)){ 
    fastq_file_path <- gsub('/$','',fastq_file_path)
    fastq_file <- vector()
    for(extension in c(".fastq.gz$", ".fastq.bz2$", ".fastq$")){
      fastq_file <- append(fastq_file, dir(fastq_file_path, extension, full.names = TRUE))
    }
  } else {fastq_file <- fastq_file_path}
  fastq_file <- normalizePath(fastq_file)
  if(length(BLAST_file_path) == 1 && dir.exists(BLAST_file_path)){ 
    BLAST_file_path <- gsub('/$','',BLAST_file_path)
    BLAST_file <- dir(BLAST_file_path, full.names = TRUE)
  } else {BLAST_file <- BLAST_file_path}
  BLAST_file <- normalizePath(BLAST_file)
  if(length(fastq_file) > 1 || length(BLAST_file) > 1){
    if(length(fastq_file) != length(BLAST_file)){
      warning('Number of submitted FASTQ and BLAST files are not equal\nonly files with corresponding names will be used.')
      fastq_file <- fastq_file[basename(sapply(strsplit(fastq_file,'\\.'),"[[", 1)) %in% basename(sapply(strsplit(BLAST_file,'\\.'),"[[", 1))]
      BLAST_file <- BLAST_file[basename(sapply(strsplit(BLAST_file,'\\.'),"[[", 1)) %in% basename(sapply(strsplit(fastq_file,'\\.'),"[[", 1))]
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

