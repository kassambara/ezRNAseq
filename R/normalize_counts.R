#' @include utilities.R
NULL
#'Normalize Counts Data Using DESeq2
#'@description Normalize counts Data using DESeq2 pipeline. Generates 3
#'  datasets: normalized count data for sequencing data; ii) rlog and iii) VST
#'  transformed data for data visualization: clustering, heatmap, PCA, ....
#'@param data an object of class data.frame, matrix, SummarizedExmeriment or
#'  DESeqDataSet containing raw counts
#'@return Return an object of class SummarizedExperiment.
#'@param thread Number of threads to be used. This depends to the available
#'  computer ressources.
#'@param result.dir a directory to save the output.
#'@param save logical value; if TRUE, the results are saved in result.dir.
#'@param size.factor.type arguments to be passed to the function
#'  DESeq2::estimateSizeFactors. Either "ratio" or "iterate".  "ratio" uses the
#'  standard median ratio method introduced in DESeq. The size factor is the
#'  median ratio of the sample over a pseudosample: for each gene, the geometric
#'  mean of all samples. "iterate" offers an alternative estimator, which can be
#'  used even when all genes contain a sample with a zero. This estimator
#'  iterates between estimating the dispersion with a design of ~1, and finding
#'  a size factor vector by numerically optimizing the likelihood of the ~1
#'  model.
#'@return returns a list containing the following components: \itemize{ \item
#'  count.norm: Normalized count for sequencing depth. The default method for
#'  estimating the size factr is "ratio".\item count.rlog: rlog-transformed data
#'  for variance stabilization. Data are in log2 scale. \item count.vst: count
#'  data after variance stabilizing transformation (vst). Data are in log2
#'  scale. \item count.fpkm: count data in fragments per kilobase per million
#'  mapped fragments.} rlog and vst data should be used for visualization such,
#'  PCA and clustering.
#'@export
normalize_counts <- function(data, result.dir = "COUNT", thread = 10,
                             save = TRUE,
                             size.factor.type = c("ratio", "iterate"))

  {
  size.factor.type <- match.arg(size.factor.type)
  BiocParallel::register( BiocParallel::MulticoreParam(workers = thread) )
  if(inherits(data, "RangedSummarizedExperiment"))
    dds <- DESeq2::DESeqDataSet(data, design = ~1)
  else if(inherits(data, "DESeqDataSet"))
    dds <- data
  else if (inherits(data, c("data.frame", "matrix"))){
    samples <- data.frame(name = colnames(data), row.names = colnames(data))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = samples,
                                          design = ~1)
  }
  else stop("Can't handle an object of class ", class(data))


  # Raw count
  # ++++++++++++++++++++++++++++
  raw.count <- as.data.frame(SummarizedExperiment::assay(dds))
  col.sums <- colSums(raw.count)
  samples.with.zero.count <- which(col.sums == 0)
  if(length(samples.with.zero.count) > 0){
    samples.names <- names(col.sums)[samples.with.zero.count]
    stop("Some samples have zero total count: ",
         paste(samples.names, collapse = ", "), ".\n",
         "Remove these samples before continuing.")
  }

  # Normalized count for sequencing depth
  # ++++++++++++++++++++++++++++
  # For bar plot
  message("- Creating normalized counts for sequencing depth...\n")
  dds <- DESeq2::estimateSizeFactors( dds, type = size.factor.type )
  count.norm <- as.data.frame(round(DESeq2::counts(dds, normalized=TRUE),2))

  # rlog-transformed data: variance stabilization
  # ++++++++++++++++++++++++++++
  # - FOR PCA, Clustering, visualization
  message("- Creating rlog data...\n")
  rld <- DESeq2::rlog(dds)
  count.rlog <- SummarizedExperiment::assay(rld)

  # Variance stabilizing transformation
  # ++++++++++++++++++++++++++++
  # - FOR PCA, Clustering, visualization
  message("- Creating variance stabilizing transformation data...\n")
  vst. <- DESeq2::varianceStabilizingTransformation(dds)
  count.vst <- SummarizedExperiment::assay(vst.)

  # FPKM
  # ++++++++++++++++++++++++++++
  message("- Creating fpkm data...\n")
  count.fpkm <- try(DESeq2::fpkm(dds))


  if(save){
    dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)
    write.table(raw.count,
                file=file.path(result.dir,  "raw.count.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
    write.table(count.norm,
                file=file.path(result.dir,  "count.normalized.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
    write.table(count.rlog,
                file=file.path(result.dir, "count.rlog.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
    write.table(count.vst,
                file=file.path(result.dir, "count.vst.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
    write.table(count.fpkm,
                file=file.path(result.dir, "count.fpkm.txt"),
                sep="\t", col.names=NA, row.names=TRUE)
  }

  res <- list(raw.count = raw.count, count.norm = count.norm,
              count.rlog = count.rlog, count.vst = count.vst, count.fpkm = count.fpkm)
  res <- structure(res, class = c("list", "normalize_counts"))
  invisible(res)
}
