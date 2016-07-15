#' Centre ranges on motif
#' 
#' @param gr A GRanges object
#' @param pwm A positional weight matrix. Should be a Matrix object with 4 rows.
#' @param bs_genome A BSgenome object.
#' @param min.score Character. Minimum alignment score percentage for the PWM.
#' @return A GRanges object. Ranges are of width=2, 
#' which are the centre of the PWM alignment.
#' @export
#' @importFrom magrittr %>%
#' @import GenomicRanges
#' @import IRanges
#' @examples
#' \dontrun{
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' gr_centre_motif(gr, ctcf.pwm, Mmusculus)
#' }
gr_centre_motif <- function(gr, pwm, bs_genome, min.score="85%"){
        
        # Check inputs
        if (class(gr) != "GRanges"){
                stop("gr is not a GRanges object!")
        }
        
        if (class(bs_genome) != "BSgenome"){
                stop("genome is not a BSgenome object!")
        }
        
        if (nrow(pwm) != 4 | class(pwm) != "Matrix"){
                stop("pwm is not correctly formated")
        }
        
        sequences <- motifRG::getSequence(gr = gr, bs_genome = bs_genome)
        
        # Function to find motif in one sequence
        find_motif_start <- function(x)
        {
                motif_starts <- Biostrings::matchPWM(pwm = pwm, subject = sequences[[x]],
                                         min.score = min.score) %>% start()
                starts <- start(gr[x]) + motif_starts
                if (length(starts) == 0){
                        out <- NULL
                } else {
                        ends <- starts + ncol(pwm)
                        out <- GenomicRanges::GRanges(seqnames = seqnames(gr[x]),
                                       ranges = IRanges(start = starts,
                                                        end = ends))
                }
                
                return(out)
        }
        
        # Apply the function across all sequences
        motif_ranges <- lapply(X = 1:length(sequences), FUN = find_motif_start)
        
        # Remove the NULLs from where no motif was detected
        is.NullOb <- function(x) is.null(motif_ranges[x]) | all(sapply(motif_ranges[x], is.null))
        no_keep <- lapply(1:length(motif_ranges), is.NullOb) %>% unlist() 
        
        # Convert to Granges object
        motif_ranges <- motif_ranges[!no_keep] %>% unlist() %>% GRangesList() %>% unlist()
        
        # Return GRanges object
        return(motif_ranges)
}


