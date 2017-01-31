#' @include utilities.R star_align.R count_reads.R
NULL
#'RNAseq Workflow
#'@description
#' RNAseq Workflow. Alignment with STAR and read counting with R/Bioconductor.
#' @inheritParams star_align
#' @inheritParams count_reads
#' @return Can create three subdirectories: \cr
#' \itemize{
#' \item SAM: containing the output of STAR alignment program
#' \item BAM: containing unsorted and sorted BAM files, and BAM index file.
#' BAM comtent is sorted by read names (*_name_sorted.bam) or by chromosome (*_sorted.bam).
#' The index files are of form *_sorted.bam.bai. *_name_sorted.bam files are used for read counting
#' using hseq-count or R. *_sorted.bam and *_sorted.bam.bai files can be used for IGV visualization.
#' \item COUNT: containing read counting results. It contains the following files:\cr
#' 1. se.RDATA: containing an object "se" which is an object of class SummarizedExperiment
#' }
#' @name rnaseq_workflow
#' @rdname rnaseq_workflow
#'@export
rnaseq_workflow <- function(data_dir = getwd(), samples.annotation = "samples.txt",
                            pairedEnd = TRUE, fastq.gz = TRUE,
                            star.index = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/StarIndex/",
                            gtf = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                            result.dir = getwd(), keep = c("name_sorted_bam"),
                            ignore.strand = FALSE, count.mode = "Union",
                            thread = 10){


  # Read samples.txt file
  samples <- read.delim(file=samples.annotation, header=TRUE, row.names=NULL)

  # Go to the directory containing the FASTQ files
  oldwd <- getwd()

  # Create result dirs
  # ++++++++++++++++++++++++++++++++++++++++
  create_dir(file.path(result.dir, "SAM"))
  create_dir(file.path( result.dir, "BAM"))
  create_dir(file.path(result.dir, "COUNT"))


  # STAR Alignement
  # ++++++++++++++++++++++++++++++++++++++++
  star_align (data_dir = data_dir, samples.annotation = samples.annotation,
              result.dir = result.dir, keep = keep,
              pairedEnd = pairedEnd, fastq.gz = fastq.gz, thread = thread,
              star.index = star.index)

  # Read counting using bioconductor
  # +++++++++++++++++++++++++++++++++++++++
  setwd(file.path(result.dir, "BAM/name_sorted"))
  se <- count_reads(bam = paste0(samples$name, "_name_sorted.bam"),
                   result.dir = file.path(result.dir, "COUNT"),
                   save = FALSE,
                   gtf = gtf, pairedEnd = pairedEnd,
                   ignore.strand = ignore.strand,
                   count.mode = count.mode, thread = thread)

  # Assign sample table to the se and save the data
  # ++++++++++++++++++++++++++++
  if("order" %in% colnames(samples)) samples <- samples[order(samples$order), ]
  if("group" %in% colnames(samples)) samples$group <- as.factor(samples$group)
  rownames(samples) <- as.vector(samples$name)
  se <- se[, as.vector(samples$name)] # same order as samples
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(samples) # add phenotypic data to se
  save(se, file = file.path(result.dir, "COUNT", "se.RDATA"))

  # Raw count & samples
  # ++++++++++++++++++++++++++++
  raw_count <- SummarizedExperiment::assay(se)
  samples <- as.data.frame(SummarizedExperiment::colData(se) )
  write.table(raw_count,
              file=file.path(result.dir, "COUNT", "raw.count.txt"),
              sep="\t", col.names=NA, row.names=TRUE)
  write.table(samples,
              file=file.path(result.dir, "COUNT", "samples.txt"),
              sep="\t", col.names=NA, row.names=TRUE)


  # Working with count data
  # ++++++++++++++++++++++++++++++++++++
  setwd(oldwd)
}
