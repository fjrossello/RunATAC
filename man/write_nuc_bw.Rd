% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_functions.R
\name{write_nuc_bw}
\alias{write_nuc_bw}
\title{Write bigwig file of nucleosome depth}
\usage{
write_nuc_bw(gr, file, scale_cpm = FALSE, nuc_span_size = c(180, 247))
}
\arguments{
\item{gr}{A GRanges object containing nucleosome positions.}

\item{file}{Path for output bigwig file. Extension should be .bw or .bigwig.}

\item{scale_cpm}{Logical. Should the output be scaled to number of 
nucleosome spanning reads.}

\item{nuc_span_size}{Numeric of length 2 specifying the insert size 
of nucleosome spanning fragnemts.}
}
\description{
Write bigwig file of nucleosome depth
}

