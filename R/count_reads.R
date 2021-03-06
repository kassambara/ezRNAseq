#' @include utilities.R
NULL
#'Read Counting
#'@description Count reads by gene and by exon
#'@param bam a vector of bam files (in the same directory) or a directory
#'  containing one or more paths to bam files. if NULL, BAM files
#'  in the current working directory are choosed.
#'@param ext BAM file extension. Default is ".bam". Used only when bam is a directory
#'@param count.mode counting mode. Read ?GenomicAlignments::summarizeOverlaps. Default
#'  value is "Union". Reads that overlap any portion of exactly one feature are
#'  counted. Reads that overlap multiple features are discarded.
#'@param result.dir a directory to save the output.
#'@param save logical value; if TRUE, the result is saved in result.dir.
#'@param by One of "gene" (for counting in gene) or "exon" (for counting in
#'  exon)
#'@param samples.annotation the path to "samples.txt" tabulation file containing sample annotations.
#'This file should contain at least 3 columns:
#'\itemize{
#'\item name: corresponding to sample name. It will be used to report the result
#'\item fastq1: the first fastq file name for paired-end sequencing
#'\item fastq2: the second fastq file name for paired-end sequencing.
#'This is optional if single-end sequencing is used
#'}
#'@param gtf gene annotation file.
#'@param pairedEnd specify if the data are paired-end sequencing data. Default is TRUE.
#'@param ignore.strand A logical indicating if strand should be considered when
#'  matching. if TRUE, allows minus strand reads to count in plus strand genes,
#'  and vice versa. Should be FALSE for stranded RNAsequencing.
#'@inheritParams fastq_nb_reads
#'@return Return an object of class SummarizedExperiment
#'@details Used Bioconductor packages: \itemize{ \item \strong{GenomicFeatures}
#'to prepare the transcript database \item \strong{Rsamtools} to import SAM/BAM
#'files in R \item \strong{GenomicAlignments} for read counting \item
#'\strong{BiocParallel}for parallel computing }
#'
#'Input files: BAM files. No need to sort the BAM files in Bioconductor version
#'> 2.12.
#'@param thread Number of threads to be used. This depends to the available computer ressources.
#'@export
count_reads <- function(bam = NULL,  ext = ".bam",
                       by = c("gene", "exon"),
                       samples.annotation = NULL,
                       result.dir = "COUNT", save = TRUE,
                       gtf = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                       pairedEnd = TRUE, ignore.strand = FALSE, count.mode = "Union", thread = 10)
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
                                             mode=count.mode, # count methods.
                                             singleEnd=singleEnd,
                                             ignore.strand = ignore.strand,
                                             fragments = fragments,
                                             inter.feature = inter.feature)
  # clearing sample names
  sple_names <- basename(bam)
  sple_names <- remove_extension(sple_names)
  sple_names <- gsub("_name_sorted", "", sple_names)
  colnames(se) <- sple_names # Rename the samples because the current name is the file name

  # Assign sample table to the se and save the data
  # ++++++++++++++++++++++++++++
  if(!is.null(samples.annotation)){
    samples <- read.delim(file=samples.annotation, header=TRUE, row.names=NULL)
    if("order" %in% colnames(samples)) samples <- samples[order(samples$order), ]
    if("group" %in% colnames(samples)) samples$group <- as.factor(samples$group)
    rownames(samples) <- as.vector(samples$name)
    se <- se[, as.vector(samples$name)] # same order as samples
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(samples) # add phenotypic data to se
  }

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
                file=file.path(result.dir, "raw.count.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
  }
  return(se)
}
