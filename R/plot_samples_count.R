#' Plot Samples Count
#' @description Plot the total count per sample
#' @param data An object of class data.frame, matrix or SummarizedExperiment.
#' @param name.size Change the sample name size
#' @inheritParams ggpubr::ggbarplot
#' @param fill fill color for the bars
#' @param ...   Add other arguments from ggbarplot to modify the plot appearance.
#' @return Return a plot ggplot.
#'
#' @export
plot_samples_count <- function(data, color = "white",
                               fill = "steelblue", ylab = "Total Count (x 10^7)",
                               xlab = "Sample", name.size = 7 , ...)

{

 # Checking the data format
  if(inherits(data, "RangedSummarizedExperiment")){
     raw_count <- SummarizedExperiment::assay(se)
  }
  else if(inherits(data, c("matrix", "data.frame"))){
    raw_count <- data
  }
  else{
    stop("wrong data format.")
  }

    res <- apply(FUN=sum, X=raw_count, MARGIN=2)

    name <- names(res)
    res_frame <- data.frame(name = factor(name, levels = name),
                            count = res/10^7)

    ggpubr::ggbarplot(res_frame, x = "name" , y = "count",
                      color = color , fill = fill,
                      xtickslab.rt = 45, font.tickslab= name.size,
                      ylab = ylab, xlab = xlab, ...)

}


