#' @include plot.R
NULL
#' Annotate using biomart
#' @description
#' \itemize{
#' \item ez_annotate(): annotate a list of gene names or ensembl ids using biomaRt
#' \item ez_biotype_stat(): Make a stat summary of gene biotype
#' }
#' @param x a vector of ids, matrix or a data frame. if x is a matrix or data frame,
#'  Then rownames is taken as ids. Ids can be Ensembl ids or gene names
#' @param dataset Dataset you want to use. See ?biomaRt::useMart
#' @name annotate
#' @rdname annotate
#' @export
ez_annotate <- function(x, host="www.ensembl.org",
                        dataset = "hsapiens_gene_ensembl")
{
  data <- NULL
  if(inherits(x, c("data.frame", "matrix"))){
    data <- x
    x <- rownames(x)
  }

  x <- unique(x)
  id_type <- .get_id_type(x)

  .load_package("biomaRt")
  ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host, dataset = dataset)
  annot <- biomaRt::getBM(mart = ensembl,
                 attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol",
                                "description", "gene_biotype",
                                "chromosome_name", "start_position", "end_position"),
                 filters = id_type,
                 values = x)
  colnames(annot) <- c("ensembl",  "entrezId", "name", "description", "biotype",
                      "chromosome", "start", "end")

  # Merge
  merge_by = "name"
  if(id_type == "ensembl_gene_id") merge_by = "ensembl"
  x <- as.data.frame(x)
  colnames(x) <- merge_by
  annot <- merge(x, annot, by.x =merge_by,
                 by.y = merge_by, all.x = TRUE)
  annot <- annot[!duplicated(as.vector(annot[, merge_by])), , drop = FALSE]
  rownames(annot) <- as.vector(annot[, merge_by])
  if(!is.null(data)){
    annot <- cbind.data.frame(annot[rownames(data), , drop = FALSE], data)
  }
  annot
}



.get_id_type <- function(id){
  n = length(id)
  id <- grep("^ENS", id)
  type <- "hgnc_symbol"
  if(length(id) == n) type <- "ensembl_gene_id"
  type
}



#' @param count_data a data frame containing raw read count
#' @param biotype the column containing biotype data of each gene
#' @param samples a vector of sample names for which you want to do the stat
#' @param compact logical for making compact statistique
#' @rdname annotate
#' @export
ez_biotype_stat <- function(count_data, biotype, samples=NULL, compact = FALSE){

  ig_gene_and_pseudo <- c('IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene', 'IG_M_gene',
                          'IG_V_gene', 'IG_Z_gene', 'IG_C_pseudogene', 'IG_J_pseudogene',
                          'IG_pseudogene', 'IG_V_pseudogene')
  protein_coding <- c('nonsense_mediated_decay', 'nontranslating_CDS', 'non_stop_decay',
                      'polymorphic_pseudogene', 'protein_coding', 'TR_C_gene', 'TR_D_gene',
                      'TR_gene', 'TR_J_gene', 'TR_V_gene')
  pseudogene <- c('disrupted_domain', 'processed_pseudogene', 'pseudogene',
                  'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene',
                  'translated_processed_pseudogene', 'translated_unprocessed_pseudogene',
                  'TR_J_pseudogene', 'TR_V_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene')

  long_non_coding_rna <- c('3prime_overlapping_ncrna', 'ambiguous_orf', 'antisense',
                           'lincRNA', 'ncrna_host', 'non_coding', 'processed_transcript',
                           'retained_intron', 'sense_intronic', 'sense_overlapping')

  short_non_coding_rna <- c('miRNA', 'miRNA_pseudogene', 'misc_RNA',
                            'misc_RNA_pseudogene', 'Mt_rRNA', 'Mt_tRNA', 'Mt_tRNA_pseudogene',
                            'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 'rRNA',
                            'rRNA_pseudogene', 'scRNA_pseudogene', 'snlRNA', 'snoRNA', 'snoRNA_pseudogene',
                            'snRNA', 'snRNA_pseudogene', 'SRP_RNA', 'tmRNA,tRNA', 'tRNA_pseudogene')


  biotype <- as.factor(count_data[, biotype])

  if(!is.null(samples)) count_data <- count_data[, samples, drop = FALSE]
  col_sum <- apply(count_data, 2, sum)
  row_sum <- apply(count_data, 1, sum)

  stat <- tapply(row_sum, biotype, sum)
  stat <- data.frame(biotype = names(stat), read_count = stat )
  stat <- stat[order(stat$biotype),]

  if(compact){

    ig_gene_and_pseudo <- sum(stat[intersect(rownames(stat), ig_gene_and_pseudo), "read_count"])
    protein_coding <- sum(stat[intersect(rownames(stat), protein_coding), "read_count"])
    pseudogene <- sum(stat[intersect(rownames(stat), pseudogene), "read_count"])
    long_non_coding_rna <- sum(stat[intersect(rownames(stat), long_non_coding_rna), "read_count"])
    short_non_coding_rna <- sum(stat[intersect(rownames(stat), short_non_coding_rna), "read_count"])
    biotype <- c("Ig_genes", "Protein_coding",
                 "Pseudogene", "Long_ncRNA", "Short_ncRNA")
    stat<-data.frame(
      biotype = factor(biotype, levels = biotype),
      read_count = c(ig_gene_and_pseudo, protein_coding, pseudogene,
                     long_non_coding_rna,  short_non_coding_rna)
    )

  }

  # total reads for all samples
  tot <- sum(col_sum)
  # plot of the % of mapped reads
  dd <- unlist(stat$read_count)
  stat$percent_mapped_reads <- dd*100/tot
  names(dd) <- substr(as.vector(stat$biotype), 1, 20)
  p <- ggbarplot(dd*100/tot, ylab="% of mapped reads")
  list(stat = stat, plot = p)
}
