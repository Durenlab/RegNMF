\name{clustering}
\alias{clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
clustering cell from result of Reg Nonnegative Matrix Factorization
}
\description{
clustering cell from result of Reg Nonnegative Matrix Factorization
}
\usage{
clustering(H)
}

\arguments{
  \item{H}{
result of RegNMF(W123$H)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ans}{Result of Rtsne}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cham
}



\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (H)
{
    UseMethod("clustering")
  }
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
