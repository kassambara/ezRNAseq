\name{ALL_bcrneg}
\alias{ALL_bcrneg}
\docType{data}
\title{B cell ALL harboring or not BCR/ABL translocation}
\description{
Microarray gene expression data from 79 individuals with B cell acute lymphoblastic leukemia (ALL). The data are RMA-normalized (log2 scale) and presented as an ExpressionSet. This data is a subset derived from the original ALL dataset availbale in the R/Bioconductor package "ALL" (https://bioconductor.org/packages/release/data/experiment/html/ALL.html).
}
\usage{data("ALL_bcrneg")}
\format{A data frame with 79 samples and 12625 probesets.}
\references{
Sabina Chiaretti, Xiaochun Li, Robert Gentleman, Antonella Vitale, Marco Vignetti, Franco Mandelli, Jerome Ritz, and Robin Foa Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. Blood, 1 April 2004, Vol. 103, No. 7.
}
\examples{
\donttest{
data(ALL_bcrneg)
dim(ALL_bcrneg)
}
}
\keyword{datasets}


