#' Read BED formatted file to GRanges object
#' 
#' @param bed_file Path Path to BED formatted file.
#' @return A GRanges object
#' @examples
#' 
#' \dontrun{
#' import_bed("myBedFile.bed")
#' }
#' @note Only imports the first 3 columns of a bed file.
read_bed <- function(bed_file){
        dat <- fread(bed_file, sep = "\t", header = FALSE)
        gr <- GenomicRanges::GRanges(seqnames = dat$V1,
                      ranges = IRanges(start = dat$V2,
                                       end = dat$V3))
        return(gr)
}


#' Read BAM file aligned fragments into GRanges object
#' 
#' @param bam_file Path to a BAM formatted file.
#' @param max_insert_size Maximum insert size allowed. Default is reccomended.
#' @return A GRanges object of aligned fragments
#' 
read_atac_frags <- function(bam_file, max_insert_size=2000){
        
        # Check inputs
        if (class(nuc_frag_size) != "numeric" | length(nuc_frag_size) !=2){
                stop("nuc_frag_size is not numeric!")
        }
        
        # Read the bam file pairs to GRanges object
        gr <- GenomicRanges::readGAlignmentPairs(file = bam_file) %>% GRanges()
        
        # Remove pairs with inserts larger than 2000 bases
        gr <- gr[width(gr) <= max_insert]
        
        return(gr)
}

#' Read BAM file Tn5 insertion positions into GRanges object
#' 
#' @param bam_file Path to a BAM formatted file.
#' @return A GRanges object of Tn5 insertion positions
#' @examples 
#' \dontrun{
#' tn <- import_atac_pos("test_100k.bam")
#' }
read_atac_pos <- function(bam_file){
        
        tn <- readGAlignments(file = bam_file) %>% GRanges()
        
        # Offset the reads to correspond to tn5 insertion site
        pos <- tn[strand(tn) == "+"] %>% 
                GenomicRanges::shift(shift=4) %>%
                resize(width = 2, fix = "start")
        strand(pos) <- "*"
        
        neg <- tn[strand(tn) == "-"] %>%
                GenomicRanges::shift(shift = -5) %>%
                resize(width = 2, fix = "start")
        strand(neg) <- "*"

        return(c(pos, neg))
}
