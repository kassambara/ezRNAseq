% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expression_set.R
\name{expression_set}
\alias{expression_set}
\alias{create_eset_tpl}
\alias{build_eset}
\title{Create an ExpressionSet}
\usage{
create_eset_tpl(exprs_file)

build_eset(exprs_file, annotation = character(), remove_tpl = FALSE)
}
\arguments{
\item{exprs_file}{path to the file containing the gene expression data}

\item{annotation}{A character describing the platform on which the samples were assayed.
This is often the name of a Bioconductor chip annotation package.}

\item{remove_tpl}{logical value. If TRUE, remove template files.}
}
\description{
Create easily an object of class ExpressionSet
using gene expression data from microarray.\cr\cr
\itemize{
\item \bold{create_eset_tpl()}: Create template files for ExpressionSet components.
\item \bold{build_eset} save ("eset.RDATA") and returns an object of class ExpressionSet (see ?Biobase::ExpressionSet)
}
}
\details{
\itemize{
 \item \bold{create_eset_tpl()}: Create template files for ExpressionSet components.
   \itemize{
   \item samples.txt: information about each sample
   \item label_description.txt: description of column names in samples.txt.
   \item experiment_description.txt: Basic description about the experiment
   }
 You must complete the content of each of these files.
 Don't forget to remove unnecessary columns.
}
}
