#' @include utilities.R summarize.R
NULL
#' Count distribution
#' @description
#' Count distribution
#' \itemize{
#' \item distribution(): plots the count distribution per sample or per group
#' }
#' @param x a numeric matrix containing raw or normalized counts; columns are samples and rows are genes
#' @param sample_cols a vector containing color for each sample.
#' @param mincount positive value; only genes with count > mincount are kept.
#' @param groupby a factor variable specifying sample groups
#' @param group_cols a vector specifying the color for each group
#' @name distribution
#' @rdname distribution
#' @export
distribution <- function(x, sample_cols = NULL, mincount = 0,
                               groupby = NULL, group_cols = NULL){


  if(!inherits(x, c("matrix", "data.frame")))
    stop("x must be a matrix or a numeric data frame")

  if(!is.null(groupby)){
    if(length(groupby) !=ncol(x))
      stop("groupby must be a factor with length equal ncol(x).")
    # Compute mean expression in each group
    x <- summarizeby(x, groupby, mean)
  }

  if(!inherits(x, "matrix")) x <- as.matrix(x)
  x <- reshape2::melt(log2(x+0.5) )
  base::colnames(x) <- c("genes", "samples", "count")
  if(mincount>0) x <- x[which(x$count >=log2(mincount)), , drop = FALSE]
  p <- ggplot2::ggplot(x, ggplot2::aes_string(x = "count", y = "..count..", color = "samples")) +
    ggplot2::geom_density()

  # Sample or group colors
  if(is.null(groupby) & !is.null(sample_cols))
    p <- p + ggplot2::scale_color_manual(values=sample_cols)
  else if (!is.null(groupby) & !is.null(group_cols))
    p <- p + ggplot2::scale_color_manual(values=group_cols)

  p <- p+
    ggplot2::labs(x = "log2(count)", y = "Frequency", title = "Count distribution")+
    ggplot2::theme_classic()

  return(p)
}

