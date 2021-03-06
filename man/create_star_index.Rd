% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_star_index.R
\name{create_star_index}
\alias{create_star_index}
\title{Generate Genome Index}
\usage{
create_star_index(fasta_file, gtf_file, thread = 4, starIndex_dir = NULL)
}
\arguments{
\item{fasta_file}{path to reference genome fasta file.
If the genome has been downloaded from \href{http://ccb.jhu.edu/software/tophat/igenomes.shtml}{iGenomes},
the path to the fasta file is of form: "~/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"}

\item{gtf_file}{path to annotation gtf file.
If iGenomes, the path can be of form: "~/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"}

\item{thread}{the number of threads to be used. Dpending on the available ressources.}

\item{starIndex_dir}{a directory to hold STAR index files. If NULL a subdirectory named "StarIndex" is created the
same directory as the fasta file.}
}
\description{
Create reference genome index for STAR software
}
\examples{
\donttest{
create_star_index(fasta_file = "~/genomes/sequence/genome.fa",
   gtf_file = "~/genomes/annotation/genes.gtf")
}

}
