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

  # Normalized count for sequencing depth
  # ++++++++++++++++++++++++++++
  # For bar plot
  message("- Writing normalized counts...\n")
  BiocParallel::register( BiocParallel::MulticoreParam(workers = thread) )
  dds <- DESeq2::DESeqDataSet(se, design = ~1)
  dds <- DESeq2::estimateSizeFactors( dds )
  count.norm <- as.data.frame(round(DESeq2::counts(dds, normalized=TRUE),2))
  write.table(count.norm,
              file=file.path(result.dir, "COUNT", "count.normalized.txt"),
              sep="\t", col.names=NA, row.names=TRUE)

  # rlog-transformed data: variance stabilization
  # ++++++++++++++++++++++++++++
  # - FOR PCA, Clustering, visualization
  message("- Writing rlog data...\n")
  rld <- DESeq2::rlog(dds)
  rld.data <- SummarizedExperiment::assay(rld)
  write.table(rld.data,
              file=file.path(result.dir, "COUNT", "count.rlog.txt"),
              sep="\t", col.names=NA, row.names=TRUE)

  # Quality Control
  # ++++++++++++++++++++++++++++
  message("- Quality Control...\n")
  nsamples <- ncol(count.norm)
  if(nsamples >= 2){
      expressed.genes <- rowSums(raw.count) > 0
      raw.count <- raw.count[expressed.genes, , drop = FALSE]
      count.norm <- count.norm[expressed.genes, , drop = FALSE]
      rld.data <- rld.data[expressed.genes, , drop = FALSE]

      # Number of mapped read counts per sample
      total.count.plot <- plot_samples_count(raw.count)
      samples.total.count <- colSums(raw.count)

      # Distribution of count per sample
      count.norm <- log2(count.norm+1)
      mcount <- tidyr::gather(count.norm, key = "samples", value = "count")
      legend. <- ifelse(nsamples > 30, "none", "bottom")
      count.dist.plot <- ggpubr::ggdensity(mcount, x = "count", color = "samples",
                                      main = "Count distribution per sample",
                                      xlab = "log2(count)", ylab = "Density",
                                      legend = legend.)


      # For heatmap and PCA
      # select genes with high sd and genes with expression > 5 in at least 10%
      expressed.genes <- genefilter::genefilter(rld.data,
                                                genefilter::pOverA(p = 0.1, 5))
      rld.data <- rld.data[expressed.genes, ]
      sds <- genefilter::rowSds(rld.data)
      top.var.genes <- head( order( sds, decreasing = TRUE ), 2000 )
      rld.data <- rld.data[top.var.genes, ]

      res.pca <- FactoMineR::PCA(t(rld.data), graph = FALSE)
      pca.plot <- factoextra::fviz_pca_ind(res.pca, repel = TRUE)

      heatmap. <- ComplexHeatmap::Heatmap( scale(t(rld.data)), name = "Exprs",
                               show_row_names = FALSE, show_column_names = TRUE,
                               show_row_dend = FALSE, show_column_dend = TRUE,
                               cluster_columns = TRUE, cluster_rows = TRUE,
                               clustering_method_columns = "complete", clustering_method_rows = "complete",
                               column_names_gp = grid::gpar(fontsize = 7),
                               column_title = "Heatmap"
                               )


      grDevices::pdf(file.path(result.dir, "COUNT", "quality.control.pdf"))
        print(total.count.plot)
        print(count.dist.plot)
        print(pca.plot)
        print(heatmap.)
      grDevices::dev.off()

      # Read Me
      sink(file.path(result.dir, "COUNT", "README.txt"))
      cat(
        "============================\n",
        "ezRNAseq R Package Workflow\n\n",
        "Data Author: ", data.author$name, " <", data.author$email, ">\n",
        "Data Analyst: ", data.analyst$name, " <", data.analyst$email, ">\n",
        "============================\n\n",
        "Date: ", Sys.Date(), "\n",
        "Alignment: STAR\n",
        "Reference Genome: ", star.index, "\n",
        "Sorting BAM Files: SAMtools\n",
        "Counting Reads: GenomicAlignments R package | Count Mode: ", count.mode, "\n\n\n",
        "Number of Mapped Reads Per Sample\n",
        "---------------------------------\n"
      )

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
