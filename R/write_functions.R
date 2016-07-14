#' Write bigwig file of nucleosome depth
#'
#' @param gr A GRanges object containing nucleosome positions.
#' @param out_file Path for output bigwig file. Extension should be .bw or .bigwig.
#' @param scale_cpm Logical. Should the output be scaled to number of 
#' nucleosome spanning reads.
#' @param nuc_span_size Numeric of length 2 specifying the insert size 
#' of nucleosome spanning fragnemts.
write_nuc_bw <- function(gr, out_file, scale_cpm=FALSE, nuc_span_size=c(180, 247)){
        
        # Get the fragment size and select only those of mononucleosome spanning size
        frag_widths <- IRanges::width(gr)
        nuc <- (frag_widths >= nuc_span_size[1] & frag_widths <= nuc_span_size[2])
        nuc_gr <- gr[nuc]
        
        # resize alignments to fragment centre to represent nucleosome position
        nuc_gr <- GenomicRanges::resize(nuc_gr, width = 2, fix = "center")
        
        # get coverage at nucleosome positions
        nuc_cov <- IRanges::coverage(x = nuc_gr)
        
        # Scale counts if specified
        if (scale_cpm == TRUE){
                nuc_cov <- (nuc_cov / length(nuc_gr)) * 1e+06
        }
        
        # write bigwig file
        rtracklayer::export.bw(object = nuc_cov, con = out_file)
}

#' Write a bigwig file of Tn5 insertion depth
#' 
#' @param gr A GRanges object containing nucleosome positions.
#' @param out_file Path for output bigwig file. Extension should be .bw or .bigwig.
#' @param scale_cpm Logical. Should the output be scaled to number of 
#' nucleosome spanning reads.
write_insertions_bw <- function(gr, out_file, scale_cpm=FALSE){
        
        # Calculate Tn insertion coverage
        cov <- IRanges::coverage(gr)
        
        # Scale counts if specified
        if (scale_cpm == TRUE){
                cov <- (cov / length(gr)) * 1e+06
        }
        
        # write Tn5 insertion bigwig
        rtracklayer::export.bw(object = cov, con = out_file)
} 


