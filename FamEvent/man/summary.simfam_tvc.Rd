\name{summary.simfam_tvc}
\alias{summary.simfam_tvc}

\title{
Summary method for \code{simfam_tvc}
}
\description{
Provides a summary of simulated data.
}
\usage{
\method{summary}{simfam_tvc}(object, digits = max(3, getOption("digits") - 3), \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object class of \code{'simfam_tvc'} generated from \code{simfam_tvc} function
}
 \item{digits}{Number of significant digits to use when printing.}
 \item{...}{Further arguments passed to or from other methods.}


}
\value{
Displays a summary of simulated data and returns the following values:
\item{num.fam}{Number of families simulated.}
\item{avg.num.affected}{Average number of affected individuals per family.}
\item{avg.num.carriers}{Average number of mutation carriers per family.}
\item{avg.family.size}{Average family size.}
\item{ave.ageonset}{Average age of onset for affected individuals.}
\item{ave.num.tvc}{Average number of TVC events per family.}
}

\author{
Yun-Hee Choi	
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{simfam_tvc}}
}
\examples{

set.seed(4321)
fam <- simfam_tvc(N.fam = 10, design = "pop", variation = "frailty", 
       base.dist = "Weibull", frailty.dist = "gamma", depend = 1, 
       add.tvc = TRUE, tvc.type = "CO", tvc.range = c(30,60), 
       tvc.parms = c(1, 0.1, 0), allelefreq = 0.02, 
       base.parms = c(0.01, 3), vbeta = c(-1.13, 2.35))

summary(fam)
\dontrun{
Study design:           pop: population-based study with affected probands 
Baseline distribution:  Weibull 
Frailty distribution:   gamma 
Number of families:                     10 
Average number of affected per family:  3.1 
Average number of carriers per family:  3.4 
Average family size:                    16.3 
Average age of onset for affected:      48.19 
Average number of TVC event per family: 4 
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }
