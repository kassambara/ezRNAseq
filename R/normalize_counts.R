#' @include utilities.R
NULL
#'Normalize Counts Data Using DESeq2
#'@description Normalize Counts Data Using DESeq2.
#'@param data an object of class data.frame, matrix, SummarizedExmeriment or
#'  DESeqDataSet containing raw counts
#'@return Return an object of class SummarizedExperiment.
#'@param thread Number of threads to be used. This depends to the available
#'  computer ressources.
#'@param result.dir a directory to save the output.
#'@param save logical value; if TRUE, the results are saved in result.dir.
#'@return returns a list containing the following components: \itemize{ \item
#'  count.norm: Normalized count for sequencing depth \item count.rlog:
#'  rlog-transformed data: variance stabilization. To be used for visualization
#'  such as PCA and clustering. }
#'@export
normalize_counts <- function(data, thread = 10, result.dir = "COUNT", save = TRUE){

  BiocParallel::register( BiocParallel::MulticoreParam(workers = thread) )
  if(inherits(data, "RangedSummarizedExperiment"))
    dds <- DESeq2::DESeqDataSet(data, design = ~1)
  else if(inherits(data, "DESeqDataSet"))
    dds <- data
  else if (inherits(data, c("data.frame", "matrix"))){
    samples <- data.frame(name = colnames(data), row.names = colnames(data))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = samples, design = ~1)
  }
  else stop("Can't handle an object of class ", class(data))


  # Raw count
  # ++++++++++++++++++++++++++++
  raw.count <- as.data.frame(SummarizedExperiment::assay(dds))

  # Normalized count for sequencing depth
  # ++++++++++++++++++++++++++++
  # For bar plot
  message("- Creating normalized counts for sequencing depth...\n")
  dds <- DESeq2::estimateSizeFactors( dds )
  count.norm <- as.data.frame(round(DESeq2::counts(dds, normalized=TRUE),2))

  # rlog-transformed data: variance stabilization
  # ++++++++++++++++++++++++++++
  # - FOR PCA, Clustering, visualization
  message("- Creating rlog data...\n")
  rld <- DESeq2::rlog(dds)
  count.rlog <- SummarizedExperiment::assay(rld)

  if(save){
    dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)
    write.table(count.norm,
                file=file.path(result.dir,  "count.normalized.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
    write.table(count.rlog,
                file=file.path(result.dir, "count.rlog.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
  }

  res <- list(raw.count = raw.count, count.norm = count.norm, count.rlog = count.rlog)
  res <- structure(res, class = c("list", "normalize_counts"))
  invisible(res)
}
