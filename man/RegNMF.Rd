\name{RgeNMF}
\alias{RgeNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Reg Nonnegative Matrix Factorization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  RgeNMF(E,
         O,
         Symbol,
         PeakName,
         Symbol_location,
         Peak_location,
         feature_cut_perc=0.01,
         corr_cut_k=100000,
         cell_num_smooth=sqrt(dim(E)[2]),
         core=8)
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{\code{matrix} or \code{sparseMatrix}. Gene by cell matrix.}
  \item{O}{\code{matrix} or \code{sparseMatrix}. Peak by cell matrix.}
  \item{Symbol}{\code{character}. Vector of genes names}

  \item{PeakName}{\code{character} Vector of peaks include its location}
  \item{Symbol_location}{\code{character} Vector of genes location}
  \item{Peak_location}{\code{character}Vector of peaks location}
  \item{feature_cut_perc}{value of filter}
  \item{corr_cut_k}{value of filter}
  \item{cell_num_smooth}{value of filter}
  \item{core}{Number of core to use}



}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{w123h$H}{H, which is result of RegNMF}
  \item{w123h$W1}{W1, which is result of RegNMF}
  \item{w123h$W2}{W2, which is result of RegNMF}
  \item{w123h$W3}{W3, which is result of RegNMF}
  \item{W123H$Reg_gene_name}{target genes of each Reg}
  \item{W123H$Reg_peak_name}{peak_name of each Reg}
  \item{W123H$Gene_name}{Reg_gene_name}
  \item{W123H$peak_name}{Reg_gene_name}
  
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Fengge Chang,\email fchang@clemson.edu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
