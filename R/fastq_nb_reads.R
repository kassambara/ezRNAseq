#' @include utilities.R
NULL
#' Count the number of reads in fastq files
#' @description Count the number of reads in fastq files.
#' @param fastq a vector or a directory containing one or more fastq files. if NULL, fastq files in the 
#' current working directory are choosed
#' @param ext fastq file extension. Default is "fastq". If you have paired-end sequencing data, you can use ext = "R1.fastq", 
#' to avoid counting twice. Note that "R1_fastq" and "R2_fastq" contain the same number of reads 
#' @param show_progress if TRUE, the progression is shown
#' @return Returns a data frame containing file name (column 1) and read count (column 2)
#' @examples
#' \donttest{
#' fastq_nb_reads(fastq, ext = "R1.fastq")
#' }
#' @name fastq_nb_reads
#' @rdname fastq_nb_reads
#' @export
fastq_nb_reads <- function(fastq = NULL, ext = "fastq", show_progress = TRUE){
  
  if(is.null(fastq)) fastq <- list_files(dir = getwd(), ext = ext, recursive = FALSE)
  else if(length(fastq) == 1 & dir.exists(fastq)) 
    fastq <- list_files(dir = fastq, ext = ext, recursive = FALSE)
  if(length(fastq) == 0) stop("No fastq files specified.")
  
  if(show_progress) cat("Start counting...\n")
  res <- NULL
  for(fq in fastq){
    bname <- basename(fq)
    if(show_progress)  cat("Counting reads in ", bname, "\n")
    read_count <- system(paste0("wc -l ", fq),  intern = TRUE)
    read_count <- strsplit(read_count, " ")[[1]]
    res <- rbind(res, c(bname, round(as.numeric(read_count[1])/4)) )
  }
  colnames(res) <- c("name", "count")
  as.data.frame(res)
}
