#' @include utilities.R
NULL
#' Generate Genome Index
#' @description
#' Create reference genome index for STAR software
#' @examples
#' \donttest{
#' create_star_index(fasta_file = "~/genomes/sequence/genome.fa",
#'    gtf_file = "~/genomes/annotation/genes.gtf")
#' }
#'
#' @param fasta_file path to reference genome fasta file.
#' If the genome has been downloaded from \href{http://ccb.jhu.edu/software/tophat/igenomes.shtml}{iGenomes},
#' the path to the fasta file is of form: "~/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
#' @param gtf_file path to annotation gtf file.
#' If iGenomes, the path can be of form: "~/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
#' @param thread the number of threads to be used. Dpending on the available ressources.
#' @param starIndex_dir a directory to hold STAR index files. If NULL a subdirectory named "StarIndex" is created the
#' same directory as the fasta file.
#'
#' @rdname create_star_index
#' @export
create_star_index <- function(fasta_file, gtf_file,  thread = 4, starIndex_dir = NULL){
  if(is.null(starIndex_dir)) starIndex_dir <- file.path(dirname(fasta_file), "StarIndex")
  create_dir(starIndex_dir)
  cmd = paste("STAR --runMode genomeGenerate --genomeDir", starIndex_dir,
              "--genomeFastaFiles", fasta_file,
              "--sjdbGTFfile", gtf_file,
              "--runThreadN", thread, sep = " ")
  system(cmd)
  cat("STAR index has been created at ", starIndex_dir, "\n\n")
}
