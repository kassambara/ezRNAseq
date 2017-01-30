#' @include utilities.R
NULL
#' Read Counting
#' @description Count reads by gene and by exon
#' @param bam a vector of bam files (in the same directory) or a directory
#'   containing one or more paths to bam files sorted by name. if NULL, BAM
#'   files in the current working directory are choosed.
#' @param ext BAM file extension. For example ext = "name_sorted.bam". Default
#'   is "name_sorted.bam". Used only when bam is a directory
#' @param mode counting mode. Read ?GenomicAlignments::summarizeOverlaps.
#'   Default value is "Union". Reads that overlap any portion of exactly one
#'   feature are counted. Reads that overlap multiple features are discarded.
#' @param result.dir a directory to save the output.
#' @param save logical value; if TRUE, the result is saved in result.dir.
#' @param by One of "gene" (for counting in gene) or "exon" (for counting in
#'   exon)
#' @inheritParams fastq_nb_reads
#' @inheritParams align
#' @return Return an object of class SummarizedExperiment
#' @details
#' Used Bioconductor packages:
#' \itemize{
#' \item \strong{GenomicFeatures} to prepare the transcript database
#' \item \strong{Rsamtools} to import SAM/BAM files in R
#' \item \strong{GenomicAlignments} for read counting
#' \item \strong{BiocParallel}for parallel computing
#' }
#'
#' Input files: BAM files. No need to sort the BAM files in Bioconductor version > 2.12
#' @export
count_reads <- function(bam = NULL,  ext = "name_sorted.bam",
                       by = c("gene", "exon"),
                       result.dir = "COUNT", save = TRUE, show_progress = TRUE,
                       gtf = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                       pairedEnd = TRUE, ignore.strand = FALSE, mode = "Union", thread = 25)
  {

  by <- match.arg(by)
  if(is.null(bam)) bam <- list_files(dir = getwd(), ext = ext, recursive = FALSE)
  else if(length(bam) == 1){
    if(dir.exists(bam)) bam <- list_files(dir = bam, ext = ext, recursive = FALSE)
  }
  if(length(bam) == 0) stop("No BAM files specified.")

  message("\n- Starting reads counting... \n")

  # 1. Prepare the gene model: Extract exons by gene
  # ================================================
  # time: ~ 10min
  message("- Preparing gene model... \n")
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
  # ================================================
  # yieldSize: control the memory. Indicates the number of reads to load at once.
  message("- Specifying bam files... \n")
  bamfiles <- Rsamtools::BamFileList(bam, yieldSize=10000000 )

  # 3. Read counting
  # ================================================
  message("- Read counting...\n")
  BiocParallel::register( BiocParallel::MulticoreParam(workers = thread) )

  singleEnd <- ifelse(pairedEnd, FALSE, TRUE)
  fragments <- ifelse(pairedEnd, TRUE, FALSE)
  se <- GenomicAlignments::summarizeOverlaps(features=features,
                                             reads=bamfiles,
                                             mode=mode, # count methods.
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
  dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)
  se_file <- file.path(result.dir, paste0("se_", by, ".RDATA"))
  if (file.exists(se_file)){
    bname <- remove_extension(basename(se_file))
    bname <- paste0(bname, "_", get_random_string(14), ".RDATA" )
    se_file2 <- file.path(dirname(se_file), bname)
    warnings("The file ", se_file, " already exists.", "The name ", se_file2, "has been used.")
    se_file <- se_file2
  }
  if(save) {
    message("- Writing count data...\n")
    save(se, file=se_file)
    raw_count <- SummarizedExperiment::assay(se)
    write.table(raw_count,
                file=file.path(result.dir, "COUNT", "raw.count.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
  }
  return(se)
}
