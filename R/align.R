#' @include utilities.R
NULL
#'RNAseq workflow using STAR and R
#'@description
#'Workflow for RNA sequencinq using STAR and R
#'@param data_dir the directory containing the data. It should contains a subdirectory called "FASTQ" which
#'contains fastq files. Default value is the working directory
#'@param samples.annotation the path to "samples.txt" tabulation file containing sample annotations.
#'This file should contain at least 3 columns:
#'\itemize{
#'\item name: corresponding to sample name. It will be used to report the result
#'\item fastq1: the first fastq file name for paired-end sequencing
#'\item fastq2: the second fastq file name for paired-end sequencing.
#'This is optional if single-end sequencing is used
#'}
#'@param star.index,gtf reference genome for STAR and annotation gtf file, respectively.
#'@param result.dir The absolute path of result directory.
#'@param keep a character vector specifying the file to be kept.
#' Allowed values include c("name_sorted_bam", "chr_sorted_bam", "sam" ).
#'@param thread Number of threads to be used. This depends to the available computer ressources.
#'@param pairedEnd specify if the data are paired-end sequencing data. Default is TRUE.
#'@param fastq.gz specify if the FASTQ files are compressed. Default is TRUE.
#'@param ignore.strand A logical indicating if strand should be considered when matching.
#'Used during read counting.
#'@param count_mode read count mode. Read ?summarizeOverlaps.
#'@details
#' The function star_align(), requires STAR and samtools programs to work.
#' Make sure that they are installed. \cr\cr
#' The workflow is as follow:\cr
#' 1. Align all FASTQ files using STAR. SAM files are generated.\cr
#' 2. Convert SAM files to BAM files. Samtools program required.\cr
#' 3. Organize BAM files (sorting and indexing). Samtools program required.\cr
#' 4. Count read using Bioconductor packages: "GenomicFeatures", "Rsamtools", "GenomicAlignments" and "BiocParallel" required.\cr
#'
#' @return Three subdirectories are created: \cr
#' \itemize{
#' \item SAM: containing the output of STAR alignment program
#' \item BAM: containing unsorted and sorted BAM files, and BAM index file.
#' BAM comtent is sorted by read names (*_name_sorted.bam) or by chromosome (*_sorted.bam).
#' The index files are of form *_sorted.bam.bai. *_name_sorted.bam files are used for read counting
#' using hseq-count or R. *_sorted.bam and *_sorted.bam.bai files can be used for IGV visualization.
#' \item COUNT: containing read counting results. It contains the following files:\cr
#' 1. se.RDATA: containing an object "se" which is an object of class SummarizedExperiment
#' }
#' @name align
#' @rdname align
#'@export
rnaseq_workflow <- function(data_dir = getwd(), samples.annotation = "samples.txt",
                       star.index = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/StarIndex/",
                       gtf = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                       result.dir = NULL,
                       keep = c("name_sorted_bam"),
                       pairedEnd = TRUE, ignore.strand = FALSE, count_mode = "Union",
                       thread = 10){


  # Read samples.txt file
  samples <- read.table(file=samples.annotation,
                        sep ="\t", header=TRUE, row.names=NULL)

  # Go to the directory containing the FASTQ files
  oldwd <- getwd()
  setwd(file.path(data_dir, "FASTQ"))

  if(is.null(result.dir)) result.dir = ".."
  else result.dir = file.path("..", result.dir)

  # Create result dirs
  # ++++++++++++++++++++++++++++++++++++++++
  create_dir(file.path(result.dir, "SAM"))
  create_dir(file.path( result.dir, "BAM"))
  create_dir(file.path(result.dir, "COUNT"))


  # STAR Alignement
  # ++++++++++++++++++++++++++++++++++++++++
  if(pairedEnd)
    cmd <- with(samples, paste("STAR --genomeDir", star.index,
                            "--readFilesIn", fastq1, fastq2,
                            "--outSAMstrandField intronMotif",
                            "--outFileNamePrefix", file.path(result.dir, "SAM", paste0(name, "_")),
                            "--runThreadN", thread, sep=" " ))
  else cmd <- with(samples, paste("STAR --genomeDir", star.index,
                            "--readFilesIn", fastq1,
                            "--outSAMstrandField intronMotif",
                            "--outFileNamePrefix", file.path(result.dir, "SAM", paste0(name, "_")),
                            "--runThreadN", thread, sep=" " ))

  for(c in cmd) system(c) # Run alignment


  # Convert SAM -> BAM file
  # ++++++++++++++++++++++++++++++++++++++++
  for(i in 1:nrow(samples)){
    sam_file <- file.path(result.dir, "SAM", paste0(samples$name[i], "_Aligned.out.sam"))
    bam_file <- file.path(result.dir, "BAM", paste0(samples$name[i], ".bam"))
    system(paste0("samtools view -bS ", sam_file, " > ", bam_file))
  }

  # Remove or keep SAM file
  # +++++++++++++++++++++++++++++++
  if(!("sam" %in% keep)){
    unlink(file.path(result.dir, "SAM"),
           recursive = TRUE, force = TRUE)
  }

  # Organize BAM files
  #++++++++++++++++++++++++++++++++++++++++++
  # - sort by name for htseq-count
  # - sort by chromosome and create index for IGV
  setwd(file.path(result.dir, "BAM"))
  for(sple in samples$name){
    bam_file <- paste0(sple, ".bam")
    system(paste0("samtools sort -n ",bam_file," ", sple,"_name_sorted"))
    if("chr_sorted_bam" %in% keep) {
      system(paste0("samtools sort ", bam_file," ", sple, "_sorted")) # by chromosome
      system(paste0("samtools index ", sple,"_sorted.bam")) # index
    }
  }


  # Read counting using bioconductor
  # +++++++++++++++++++++++++++++++++++++++
  se <- count_reads(bam = paste0(samples$name, "_name_sorted.bam"),
                   result.dir = file.path(result.dir, "COUNT"),
                   save = FALSE,
                   show_progress = TRUE,
                   gtf = gtf, pairedEnd = pairedEnd,
                   ignore.strand = ignore.strand,
                   mode = count_mode, thread = thread)

  # Assign sample table to the se and save the data
  # ++++++++++++++++++++++++++++
  if("order" %in% colnames(samples)) samples <- samples[order(samples$order), ]
  if("group" %in% colnames(samples)) samples$group <- as.factor(samples$group)
  rownames(samples) <- as.vector(samples$name)
  se <- se[, as.vector(samples$name)] # same order as samples
  GenomicRanges::colData(se) <- S4Vectors::DataFrame(samples) # add phenotypic data to se
  save(se, file = file.path(result.dir, "COUNT", "se.RDATA"))

  # Raw count & samples
  # ++++++++++++++++++++++++++++
  raw_count <- SummarizedExperiment::assay(se)
  samples <- as.data.frame(SummarizedExperiment::colData(se) )
  write.table(raw_count,
              file=file.path(result.dir, "COUNT", "raw.count.txt"),
              sep="\t", col.names=NA, row.names=TRUE)
  write.table(samples,
              file=file.path(result.dir, "COUNT", "samples.txt"),
              sep="\t", col.names=NA, row.names=TRUE)


  # Working with count data
  # ++++++++++++++++++++++++++++++++++++
  setwd(oldwd)
}

#' @export
#' @describeIn align Aligning RNAseq data using STAR software
star_align <- function(data_dir = getwd(), samples.annotation = "samples.txt",
                       star.index = "/eqmoreaux/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/StarIndex/",
                       result.dir = NULL,
                       keep = c("name_sorted_bam"),
                       pairedEnd = TRUE, fastq.gz = TRUE, thread = 10)
{

  # Read samples.txt file
  samples <- read.delim(file=samples.annotation, header = TRUE, row.names = NULL)

  # Go to the directory containing the FASTQ files
  oldwd <- getwd()
  setwd(file.path(data_dir, "FASTQ"))

  if(is.null(result.dir)) result.dir = ".."

  # Create result dirs
  # ++++++++++++++++++++++++++++++++++++++++
  create_dir(file.path(result.dir, "SAM"))
  create_dir(file.path( result.dir, "BAM"))


  # STAR Alignement
  # ++++++++++++++++++++++++++++++++++++++++
  message("\nStarting Alignment...\n",
          "======================\n")

  # Option to read compressed FASTQ files
  if(fastq.gz) read_file_command <- "--readFilesCommand gunzip -c"
  else read_file_command <- ""

  if(pairedEnd)
    cmd <- with(samples, paste("STAR --genomeDir", star.index,
                               read_file_command,
                               "--readFilesIn", fastq1, fastq2,
                               "--outSAMstrandField intronMotif",
                               "--outFileNamePrefix", file.path(result.dir, "SAM", paste0(name, "_")),
                               "--runThreadN", thread, sep=" " ))
  else cmd <- with(samples, paste("STAR --genomeDir", star.index,
                                  read_file_command,
                                  "--readFilesIn", fastq1,
                                  "--outSAMstrandField intronMotif",
                                  "--outFileNamePrefix", file.path(result.dir, "SAM", paste0(name, "_")),
                                  "--runThreadN", thread, sep=" " ))

  for(c in cmd) {
    cat(c, "\n--------------------------------------\n")
    system(c) # Run alignment
  }
  message("\n####Alignment Finished###\n")



  # Convert SAM -> BAM file
  # ++++++++++++++++++++++++++++++++++++++++
  message("\nConverting SAM to BAM Files...\n",
          "======================\n")
  sam_file <- file.path(result.dir, "SAM/{}_Aligned.out.sam")
  bam_file <- file.path(result.dir, "BAM/{}.bam")
  cmd <- paste(paste(samples$name), "|", "parallel -j", thread,
             "samtools view -bs", sam_file, ">", bam_file,
              sep = " ")
  system(cmd)

  # Remove or keep SAM file
  # +++++++++++++++++++++++++++++++
  if(!("sam" %in% keep)){
    unlink(file.path(result.dir, "SAM"),
           recursive = TRUE, force = TRUE)
  }

}
