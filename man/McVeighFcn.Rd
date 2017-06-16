\name{McVeighFcn}
\alias{McVeighFcn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{McVeighFcn
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
McVeighFcn is used to process minute-level accelerometer data from devices worn over 24 hours with the McVeigh algorithm to identify instances of in-bed versus sleep time.  The McVeigh specification for determining non-wear time can be used or a non-wear vector calculated via some other methodology can be input.
}
\usage{
McVeighFcn(timevec, cpmvec, nwvec= NULL, rollwin = 30,
slthresh = 88, slperiod = 180, wkperiod = 10,
prslthresh1 = 89, prslthresh2 = 50, prslcount = 4,
prwkthresh1 = 91, prwkthresh2 = 200, prwkcount = 3,
spurthresh = 2.5*60, nwperiod = 90, nwctsgt0 = 3,
nwmaxct = 50,    vhighpct = 0.9, highpct = 0.8,
modpct = 0.5,fairtime = 14*60,  validday = 10*60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{timevec}{vector detailing the time of each observation.  Must be minute-level data.  The specified format is “POSIX”.  When \code{timevec} is a “POSIX” variable, an error will be returned if successive measurements are not separated by one minute.  If \code{timevec} is not a “POSIX” variable, rows will be assumed to spaced one-minute apart and the calculation of \code{validday} and \code{validinbed} variables will not be possible.}
\item{cpmvec}{vector summarizing activity level in cpm at each observation point. Must be same length as \code{timevec}.}
\item{nwvec}{optional vector specifying non-wear observations.  If left as NULL, nonwear observations will be calculated via McVeigh’s algorithm.}
\item{rollwin}{number of minutes on either side of an observation to be used to calculate moving/rolling average.  If the entire size of \code{rollwin} is not available (e.g., observation 2, with \code{rollwin} = 5), the number number of measurements available are used (e.g, observations 1:8 in the example above).}
\item{slthresh}{threshold defining low versus high activity.}
\item{slperiod}{consecutive time (in minutes) of rolling average low activity require to define approximate into-bed times.}
\item{wkperiod}{consecutive time (in minutes) of rolling average high activity require to define approximate out-of-bed times.}
\item{prslthresh1}{threshold for individual cpm measurements when calculating precise into-bed times.}
\item{prslthresh2}{along with \code{prslcount}, used to assess 10-minute intervals to calculate precise into-bed times.}
\item{prslcount}{threshold for number of cpm observations below \code{prslthresh2} when calculating precise into-bed times.}
\item{prwkthresh1}{threshold for individual cpm measurements when calculating precise out-of-bed times.}
\item{prwkthresh2}{Along with \code{prwkcount}, used to assess 10-minute intervals to calculate precise out-of-bed times.}
\item{prwkcount}{threshold for number of cpm observations below \code{prwkthresh2} when calculating precise out-of-bed times.}
\item{spurthresh}{time in minutes between consecutive in-bed periods used to classify spurious data.  When the interval between 2 in-bed periods  < \code{spurthresh}, a new in-bed period corresponding to this interval is created.}
\item{nwperiod}{length of interval in minutes considered for non-wear calculations.}
\item{nwctsgt0}{number non-zero measurements within \code{nwperiod} allowed for non-wear classification.}
\item{nwmaxct}{threshold for maximum value allowed when classifying non-wear intervals.}
\item{vhighpct}{in-bed period with a proportion ≥ \code{vhighpct} of non-wear measurements is considered a non-wear in-bed period.}
\item{highpct}{in-bed period with a proportion ≥ \code{highpct} of non-wear measurements and remaining proportion ≤ \code{slperiod} is considered a non-wear in-bed period.}
\item{modpct}{in-bed period with a proportion ≤ \code{modpct} of non-wear measurements that is less than ≤ \code{fairtime} is considered a wear in-bed period.}
\item{fairtime}{maximum time allowed for wear in-bed period.}
\item{validday}{number of minutes of waking wear required for a day to be considered valid.}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
McVeighFcn returns a data frame with the following variables:

\item{time}{the input \code{timevec} vector.}
\item{cpm}{the input \code{cpmvec} vector.}
\item{nw}{either the input \code{nwvec} vector or the calculated non-wear vector.}
\item{inbed}{binary variable indicating whether observation is an in-bed observation.}
\item{wakingwear}{binary variable indicating whether observation is a waking wear observation.}
\item{validday}{binary variable indicating whether observation occurred on a valid day.  If \code{timevec} is not of type “POSIX”, \code{NA} will be returned for this variable. }
\item{validinbed}{binary variable indicating whether observation is an inbed observation from a valid day.  If \code{timevec} is not of type “POSIX”, NA will be returned for this variable. }

}
\references{
McVeigh, J. A., Winkler, E. A., Healy, G. N., Slater, J., Eastwood, P. R., & Straker, L. M. (2016). Validity of an automated algorithm to identify waking and in-bed wear time in hip-worn accelerometer data collected with a 24 h wear protocol in young adults. \emph{Physiological Measurement}, 37(10), 1636.
}
\author{
Vincent Berarardi,Assistant Professor, Chapman University, Department of Psychology
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
