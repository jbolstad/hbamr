\name{LC1980}
\alias{LC1980}

\title{
  1980 Liberal-Conservative Scales
}

\description{
  Liberal-Conservative 7-point scales from the 1980 National Election Study. Includes
  (in order) self-placement, and rankings of Carter, Reagan, Kennedy, Anderson,
  Republican party, Democratic Party. Stored as a matrix of integers. The numbers 0,
  8, and 9 are considered to be missing values.
}

\usage{
data(LC1980)
}

\value{

  The data is formatted as an integer matrix with the following elements.

  \item{LC1980}{ matrix, containing reported placements of various stimuli on a
	7 point Liberal-Conservative scale:
    \itemize{
    \item{\code{Self}}{ Self-placement on 7 point scale.}
    \item{\code{Carter}}{ Placement of Carter on 7 point scale.}
    \item{\code{Reagan}}{ Placement of Reagan on 7 point scale.}
    \item{\code{Kennedy}}{ Placement of Kennedy on 7 point scale.}
    \item{\code{Anderson}}{ Placement of Anderson on 7 point scale.}
    \item{\code{Republicans}}{ Placement of Republicans on 7 point scale.}
    \item{\code{Democrats}}{ Placement of Democrats on 7 point scale.}
    }
  }
}

\source{
  This dataset was originally part of the `basicspace` package.

  American National Election Studies. \url{http://www.electionstudies.org/}

  Also availble from Keith Poole. \url{ http://voteview.com }
}

\author{
Keith Poole \email{ktpoole@uga.edu}

Howard Rosenthal \email{hr31@nyu.edu}

Jeffrey Lewis \email{jblewis@ucla.edu}

James Lo \email{lojames@usc.edu}

Royce Carroll \email{rcarroll@rice.edu}
}

\examples{
### Loads Liberal-Conservative scales from the 1980 NES.
data(LC1980)
dat <- LC1980
dat[dat == 0 | dat == 8 | dat == 9] <- NA
self <- dat[, 1]
stimuli <- dat[, -1]
}

\keyword{datasets}
