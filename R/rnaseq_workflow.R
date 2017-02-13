#' @include utilities.R star_align.R count_reads.R
NULL
#'RNAseq Workflow
#'@description
#' RNAseq Workflow. Alignment with STAR and read counting with R/Bioconductor.
#' @inheritParams star_align
#' @inheritParams count_reads
#' @param data.analyst a list containing the name and the email of the data analyst.
#' @param data.author a list containing the name and the email of the data author.
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
                            thread = 10,
                            data.analyst = list(name = "", email = ""),
                            data.author = list(name = "", email = "")
                            )
  {


  # Read samples.txt file
  samples <- read.delim(file=samples.annotation, header=TRUE, row.names=NULL)

  # Go to the directory containing the FASTQ files
  oldwd <- getwd()

  # Create result dirs
  # ++++++++++++++++++++++++++++++++++++++++
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
  message("- Writing raw counts...\n")
  raw.count <- as.data.frame(SummarizedExperiment::assay(se))
  samples <- as.data.frame(SummarizedExperiment::colData(se) )
  write.table(raw.count,
              file=file.path(result.dir, "COUNT", "raw.count.txt"),
              sep="\t", col.names=NA, row.names=TRUE)
  write.table(samples,
              file=file.path(result.dir, "COUNT", "samples.txt"),
              sep="\t", col.names=NA, row.names=TRUE)

  # - Normalizing count for sequencing depth: can be used for bar plots
  # - rlog-transformed data: variance stabilization. FOR PCA, Clustering, visualization
  # ++++++++++++++++++++++++++++
  count.data <- try(normalize_counts(se, thread = thread,
                                  result.dir = file.path(result.dir, "COUNT"), save = TRUE))

  # Quality Control
  # ++++++++++++++++++++++++++++
  if(!inherits(count.data, "try-error"))
    check_count_data(count.data, result.dir = file.path(result.dir, "COUNT"))

  nsamples <- ncol(raw.count)
  if(nsamples >= 2){

      # Read Me
      sink(file.path(result.dir, "COUNT", "README.txt"))
      cat(
        "============================\n",
        "ezRNAseq R Package Workflow\n\n",
        "Data Author: ", data.author$name, " <", data.author$email, ">\n",
        "Data Analyst: ", data.analyst$name, " <", data.analyst$email, ">\n",
        "============================\n\n",
        "Date: ", as.character(Sys.Date()), "\n",
        "Alignment: STAR\n",
        "Reference Genome: ", star.index, "\n",
        "Sorting BAM Files: SAMtools\n",
        "Counting Reads: GenomicAlignments R package | Count Mode: ", count.mode, "\n\n\n",
        "Number of Mapped Reads Per Sample\n",
        "---------------------------------\n"
      )
      # Number of mapped read counts per sample
      samples.total.count <- colSums(raw.count)

      df <- data.frame(
        samples = names(samples.total.count),
        total_mapped_reads = samples.total.count
      )
      rownames(df) <- 1:nrow(df)
      print(df)
      sink()
  }


  # Working with count data
  # ++++++++++++++++++++++++++++++++++++
  setwd(oldwd)
}
