#' Read a BAM file and get tn5 positions and create a BED file centered on tn5 insertion sites
#' @param bam_file Path to a BAM formatted file.
#' @param bed_file Path for output BED file
#' @param width The size of the interval to be centered on the tn5 insertion site
#' @return No value is returned. 
#' \dontrun{
#' bam_to_offset_bed(bam_file="test_100k.bam", bed_file="test_100k_offset.bed")
#' }
#' @export
bam_to_offset_bed <- function(bam_file, bed_file, width=50){
        gr <- read_atac_pos(bam_file)
        gr <- GenomicRanges::resize(gr, width = width, fix = "center")
        write_gr_bed(gr = gr, file = bed_file)
}