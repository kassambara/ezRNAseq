#' @include utilities.R
NULL
#' Download Reference Genome from iGenome
#' @description
#' Download and install illumina \href{http://ccb.jhu.edu/software/tophat/igenomes.shtml}{iGenomes}.
#'
#' @param version the version of the genome of interest (e.g., "hg19"). Allowed values are "hg18", "hg19" and "GRCh37".
#' For installing other iGenome versions, specify the path
#' @param path the path to the genome of interest. Not required if version = "hg18", "hg19" or "GRCh37".
#' @param result_dir a directory to hold the genome. Default is "~/genomes".
#' @return Create a result directory containing the reference igenome.
#' Default directory is "~/genomes" which will be created in your home directory.
#' @examples
#' \donttest{
#' install_igenome(version = "grch37", result_dir = "~/genomes")
#' }
#'
#' @rdname install_igenome
#' @export
install_igenome <- function(version = NULL, path = NULL, result_dir = "~/genomes"){

  if(is.null(version) & is.null(path))
    stop("Specify an iGenome version or path.")
  else if(!is.null(version)){
    .check_igenome_version(version)
    path <- .build_igenome_path(version)
  }

  cat("Installing iGenome\n - version: ", version, "\n",
      "- Path: ", path, "\n\n")

  oldwd <- getwd()
  create_dir(result_dir)
  setwd(result_dir)
  system(paste0("wget ", path)) # download
  system(paste0("tar zxvf ", basename(path) )) # uncompress
  system(paste0("rm ", basename(path) ))

  setwd(oldwd)
}


.check_igenome_version <- function(version){
  version <- tolower(version)
  if(!(version %in% c("hg18", "hg19", "grch37")))
    stop("Can't handle a genome version ", version,
         ". Please specify an iGenome path choosed from http://ccb.jhu.edu/software/tophat/igenomes.shtml")
}

.build_igenome_path <- function(version){

  ftp_igenome <- "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com"
  version <- tolower(version)
  if(version %in% c("hg18", "hg19")){
    path <- paste0(ftp_igenome, "/Homo_sapiens/UCSC/",
                  version, "/Homo_sapiens_UCSC_",
                  version, ".tar.gz" )
  }
  else if(version == "grch37"){
    path <- paste0(ftp_igenome, "/Homo_sapiens/Ensembl/",
                   "GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz")
  }
  else stop("Can't handle a genome version ", version)
}
