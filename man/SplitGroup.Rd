\name{SplitGroup}
\alias{SplitGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SplitGroup
}
\description{
Allot cells(barcords) and regulations(peak - gene) to each clusters
}
\usage{

}

\arguments{
  \item{foldername}{output  folder}
  \item{barcord}{barcode of cell}
  \item{W3}{W3 in W123H}
  \item{H}{H in W123H}
  \item{Reg_symbol_name}{Gene_name in W123H}
  \item{Reg_peak_name}{peak_name in W123H}
  \item{cluster}{}
}
\details{
About argument barcord and cluster, should be matched.
For example, if there are barcode A, B, C, belong cluster 1, 2, 1, 
Then argument `barcode` and `cluster` should be \code{Vector} (A,B,C) and \code{Vector}(1,2,1)
}
\value{
  \item{barcordFileName}{A tsv file contents barcodes with their cluster information}
  \item{RegFolderName}{A folder contents regulations in each cluster as tsv file, will be filtered by peak of each cluster}
}


\author{
Cham
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

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
