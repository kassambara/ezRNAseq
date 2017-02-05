#'Create Samples Annotation File
#'@description Create samples annotation file according fastq files contained in a FASTQ directory.
#'@param fastq_dir path to the directory containing the fastq files.
#'@param fastq1_pattern character string in the fastq1 files
#'@param fastq2_pattern character string in the fastq2 files (for paired-end sequencing).
#'@return Create a samples.txt file in the working directory.
#'@export
create_samples_file <- function(fastq_dir = "FASTQ",
                                fastq1_pattern = "_R1.fastq.gz",
                                fastq2_pattern = "_R2.fastq.gz")
{
  fastq1_files <- list.files(fastq_dir, pattern = fastq1_pattern)
  fastq2_files <- list.files(fastq_dir, pattern = fastq2_pattern)
  pairedend <- TRUE

  if(length(fastq1_files) == 0)
    stop("Can't find any fastq1 files matching the provided pattern")

  if(length(fastq2_files) == 0){
    warning("Can't find any fastq2 files matching the provided pattern. ",
            "I guess this is a single-end sequencing data.")
    pairedend <- FALSE
  }

  name1 <- gsub(fastq1_pattern, "", fastq1_files)
  samples.data <- data.frame(
    name = name1,
    fastq1 = fastq1_files
  )


  if(pairedend){
    name2 <- gsub(fastq1_pattern, "", fastq1_files)
    if(!all(name1==name2))
      stop("FASTQ1 sample names are different from FASTQ2 sample names.")
    samples.data$fastq2 = fastq2_files
  }

  write.table(samples.data, file = "samples.txt",
              sep="\t", row.names=FALSE)

  message("Sample annotation file saved at: ", file.path(getwd(), "samples.txt\n"))
}
