#' @include utilities.R
NULL
#' Create an ExpressionSet
#'@description
#' Create easily an object of class ExpressionSet
#' using gene expression data from microarray.\cr\cr
#' \itemize{
#' \item \bold{create_eset_tpl()}: Create template files for ExpressionSet components.
#' \item \bold{build_eset} save ("eset.RDATA") and returns an object of class ExpressionSet (see ?Biobase::ExpressionSet)
#' }
#' @param exprs_file path to the file containing the gene expression data
#' @param annotation A character describing the platform on which the samples were assayed.
#' This is often the name of a Bioconductor chip annotation package.
#' @details
#' \itemize{
#'  \item \bold{create_eset_tpl()}: Create template files for ExpressionSet components.
#'    \itemize{
#'    \item samples.txt: information about each sample
#'    \item label_description.txt: description of column names in samples.txt.
#'    \item experiment_description.txt: Basic description about the experiment
#'    }
#'  You must complete the content of each of these files.
#'  Don't forget to remove unnecessary columns.
#' }
#' @name expression_set
#' @rdname expression_set
#'@export
create_eset_tpl <- function(exprs_file){
  data_dir <- dirname(exprs_file)
  exprs <- read.table(exprs_file, sep="\t", header = TRUE,
                      row.names = 1, as.is = TRUE)

  # Create expression matrix
  write.table(exprs, file = file.path(data_dir, "exprs.txt"), sep ="\t",
              row.names = TRUE, col.names = NA, quote = FALSE)

  # Create phenotypic data: sample annotation
  nsples <- ncol(exprs)
  sples <- data.frame(
    name = rep("---", nsples), title = rep("---", nsples), sample_order = 1:nsples,
    group = rep("---", nsples), group_title = rep("---", nsples), group_order = rep("---", nsples),
    etc... = rep("---", nsples), row.names = colnames(exprs))
  write.table(sples, file = file.path(data_dir, "samples.txt"), sep ="\t",
              row.names = TRUE, col.names = NA, quote = FALSE)

  # Experiment description
  label_desc <- data.frame( labelDescription = c("sample short name",
                         "long name of sample", "sample order",
                         "sample group name", "group title", "group order"),
    row.names = c("name", "title", "sample_order", "group", "group_title", "group_order"))
  write.table(label_desc, file = file.path(data_dir, "label_description.txt"), sep ="\t",
              row.names = TRUE, col.names = NA, quote = FALSE)

  # Experiment description
  experiment_desc <- data.frame(
    row.names = c("name", "lab", "contact", "title", "abstract", "url", "PMIDs", "etc..."),
    value = c("Pierre Fermat", "Francis Galton Lab", "pfermat@lab.not.exist", "Smoking-Cancer Experiment",
              "An example ExpressionSet", "www.lab.not.exist", "pubmed id if published",  "etc ..."))
  write.table(experiment_desc, file = file.path(data_dir, "experiment_description.txt"), sep ="\t",
              row.names = TRUE, col.names = NA, quote = FALSE)
  message("The template files (",
          "exprs.txt, samples.txt, label_description.txt and experiment_description.txt) ",
          "has been created in the folder: ", data_dir)
}

#' @rdname expression_set
#' @param remove_tpl logical value. If TRUE, remove template files.
#'@export
build_eset <- function(exprs_file, annotation = character(), remove_tpl = FALSE){
  data_dir <- dirname(exprs_file)
  # Import expression file as a matrix
  exprs <- as.matrix(read.table(exprs_file, header = TRUE,
                                sep = "\t", row.names = 1, as.is = TRUE))
  # Phenotypic data
  # ++++++++++++++++++
  pData <- read.table(file.path(data_dir, "samples.txt"), row.names = 1, header = TRUE,
                      sep = "\t")
  if(!all(rownames(pData) == colnames(exprs)))
    stop("The number and name of rows in 'samples.txt' ",
         "must match the number and names of columns in expression file.")

  # Label description
  label_desc <- read.table(file.path(data_dir, "label_description.txt"),
                           row.names = 1, header = TRUE, sep = "\t")
  if(!all(rownames(label_desc) == colnames(pData)))
    stop("The row names of 'label_description.txt' must be identical ",
         "to the column names of 'samples.txt'.")
  # Create phenoData
  phenoData <- new("AnnotatedDataFrame", data = pData, varMetadata = label_desc)

  # Experiment description
  # ++++++++++++++
  exp_desc <- as.matrix(read.table(file.path(data_dir, "experiment_description.txt"),
                           row.names = 1, header = TRUE, sep = "\t"))
  d <- setdiff(rownames(exp_desc), c("name", "lab", "contact", "title", "abstract", "url"))
  others = list()
  if(length(d)> 0){
    for(i in 1:length(d)) others[[i]] <- exp_desc[d[i], 1]
    names(others) <- d
  }
  experimentData <- new("MIAME",
                        name = exp_desc["name", 1],
                        lab = exp_desc["lab", 1],
                        contact=exp_desc["contact", 1],
                        title=exp_desc["title", 1],
                        abstract=exp_desc['abstract', 1],
                        url=exp_desc["url", 1],
                        other=others)

  # Assembling
  # +++++++++++++++++
  eset <- Biobase::ExpressionSet(assayData=exprs,
                              phenoData=phenoData,
                              experimentData=experimentData,
                              annotation= annotation)

  if(remove_tpl) file.remove(file.path(data_dir, c("samples.txt", "label_description.txt",
                                                   "experiment_description.txt")))

  save(eset, file = file.path(data_dir, "eset.RDATA"))
  eset
}
