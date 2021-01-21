#' Find a gene or a nearest gene for genomic positions or genomic ranges
#'
#'
#' @param chr a vector (or a value) of chromosome numbers for genomic ranges.
#' @param start a vector (or a value) specifying start positions for genomic ranges.
#' @param end a vector (or a value) specifying end positions for genomic ranges. If NULL, start will be used. Default is NULL.
#' @param GRCh a value (37 or 38) indicating which reference genome will be used. GRCh 37 and GRCh38 only available. Default is GRCh=37.
#'
#' @examples
#' install.packages("BiocManager")
#' BiocManager::install(c("ensembldb","EnsDb.Hsapiens.v75","EnsDb.Hsapiens.v86"))
#' chr = rep(21,5)
#' start = c(36470216,36470077,36471055,36470098,36470139)
#' end = start + 10
#' locate_gene(chr=chr,start=start,end=end)
#'
#' @author Hyo Young Choi
#' @import ensembldb EnsDb.Hsapiens.v75 EnsDb.Hsapiens.v86
locate_gene = function(chr,start,end=NULL,GRCh=37) {
  require(ensembldb)

  if (GRCh==37) {
    require(EnsDb.Hsapiens.v75) # GRCh37
    edb = EnsDb.Hsapiens.v75  # abbreviate
  } else if (GRCh==38) {
    require(EnsDb.Hsapiens.v86) # GRCh38
    edb = EnsDb.Hsapiens.v86  # abbreviate
  } else {
    stop(paste0("Reference genome GRCh ",GRCh," not available. Please use 37 or 38."))
  }
  chr = as.numeric(chr)
  start = as.numeric(start)
  if (is.null(end)) {
    end = start
  } else {
    end = as.numeric(end)
  }
  # Select protein coding
  Tx = transcripts(edb,
                   columns = c(listColumns(edb , "tx"), "seq_name", "gene_name", "gene_seq_start", "gene_seq_end"),
                   filter = TxBiotypeFilter("protein_coding"),
                   return.type = "DataFrame")

  # Find genes containing the region
  gr = GRanges(chr, ranges = IRanges(start,end))
  grf = GRangesFilter(gr, type = "any")
  overlap = genes(edb, filter = ~ tx_biotype %in% c("protein_coding") & grf)

  # Find the nearest gene (k=1)
  # Make Tx GRanges object so that it can be used in "nearest" function.
  Tx.gr = GRanges(seqnames=Tx$seq_name,ranges=IRanges(start=Tx$gene_seq_start,end=Tx$gene_seq_end))
  # Also, make another GRange object containing the region of interest.
  grf2 = GRangesFilter(Tx.gr[nearest(gr,Tx.gr)])
  near.overlap = genes(edb, filter = ~ tx_biotype %in% c("protein_coding") & grf2)

  return(list(Gene=overlap$symbol,Nearest.Gene=near.overlap$symbol[which(! near.overlap$symbol %in% overlap$symbol)]))
}
