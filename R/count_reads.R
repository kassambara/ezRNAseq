#' @include utilities.R
NULL
#' Count reads per genes using Bioconductor package
#' @param bam a vector of bam files (in the same directory)
#'  or a directory containing one or more paths to bam files sorted by name.
#' if NULL, BAM files in the current working directory are choosed.
#' @param ext BAM file extension. For example ext = "name_sorted.bam".
#' Default is "name_sorted.bam". Used only when bam is a directory
#' @param mode counting mode. Read ?summarizeOverlaps
#' @param result_dir a directory to save the output.
#' @param save logical value; if TRUE, the result is saved in result_dir.
#' @param by One of "gene" (for counting in gene) or "exon" (for counting in exon)
#' @inheritParams fastq_nb_reads
#' @inheritParams align
#' @return
#' Return an object of class SummarizedExperiment
#' @details required bioconductor packages: "GenomicFeatures", "Rsamtools", "GenomicAlignments" and "BiocParallel".
#' @export
count_reads <- function(bam = NULL,  ext = "name_sorted.bam",
                       by = c("gene", "exon"),
                       result_dir = "COUNT", save = TRUE, show_progress = TRUE,
                       gtf = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                       pairedEnd = TRUE, ignore.strand = FALSE, mode = "Union", thread = 25)
  {

  by <- match.arg(by)
  if(is.null(bam)) bam <- list_files(dir = getwd(), ext = ext, recursive = FALSE)
  else if(length(bam) == 1){
    if(dir.exists(bam)) bam <- list_files(dir = bam, ext = ext, recursive = FALSE)
  }
  if(length(bam) == 0) stop("No BAM files specified.")

  if(show_progress) cat("Checking packages for read counting... \n")
  # Check packages and install missing ones
  pkgs <- c("GenomicFeatures", "Rsamtools", "GenomicAlignments", "BiocParallel")
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  install_pkgs(pkgs_miss, check = FALSE)

  if(show_progress) cat("Preparing gene model... \n")
  # 1. Prepare the gene model: Extract exons by gene
  requireNamespace("GenomicFeatures", quietly = TRUE)
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format="gtf" )
  if(by == "gene") {
    features <- GenomicFeatures::exonsBy(txdb, by="gene" )
    inter.feature <- TRUE
  }
  else if(by == "exon") {
    features <- GenomicFeatures::disjointExons(txdb)
    inter.feature <- FALSE
  }

  # 2. Specify bam files
  requireNamespace("Rsamtools", quietly = TRUE)
  bamfiles <- Rsamtools::BamFileList(bam, yieldSize=10000000 )


  # 3. Read counting
  if(show_progress) cat("Read counting...\n")
  requireNamespace("GenomicAlignments", quietly = TRUE)
  requireNamespace("BiocParallel", quietly = TRUE)
  BiocParallel::register( BiocParallel::MulticoreParam(workers = thread) )
  singleEnd = TRUE
  fragments = FALSE
  if(pairedEnd) {
    singleEnd = FALSE
    fragments = TRUE
  }
  se <- GenomicAlignments::summarizeOverlaps(features=features,
                                             reads=bamfiles,
                                             mode=mode,
                                             singleEnd=singleEnd,
                                             ignore.strand = ignore.strand,
                                             fragments = fragments,
                                             inter.feature = inter.feature)
  # clearing sample names
  sple_names <- basename(bam)
  sple_names <- remove_extension(sple_names)
  sple_names <- gsub("_name_sorted", "", sple_names)
  colnames(se) <- sple_names # Rename the samples because the current name is the file name
  # save se file
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
  se_file <- file.path(result_dir, paste0("se_", by, ".RDATA"))
  if (file.exists(se_file)){
    bname <- remove_extension(basename(se_file))
    bname <- paste0(bname, "_", get_random_string(14), ".RDATA" )
    se_file2 <- file.path(dirname(se_file), bname)
    warnings("The file ", se_file, " already exists.", "The name ", se_file2, "has been used.")
    se_file <- se_file2
  }
  if(save) save(se, file=se_file)
  return(se)
}
