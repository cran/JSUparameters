\name{JSUparameters}
\alias{JSUparameters}
\title{Find the Best-Fitting JohnsonSU Distribution to a dataset.
}
\description{
This function finds the best-fitting JohnsonSU distribution for a given dataset through least squares optimisation.
}
\usage{
JSUparameters(dat)
}
\arguments{
  \item{dat}{
The data you would like to fit a JohnsonSU distribution to.
}
}
\details{
This function contains many sub-funtcions.
1. OptimiseGivenZeta(): This function will find the best-fitting JohnsonSU distribution for the given dataset, for a given value of zeta (typically found using a GSS).
2. golden(): Carries out a Golden Section Search (GSS) to find the minimum of a function within a given interval.
3. calculate_ssq(): Calculates the sum of squares (SSQ) for a given value of zeta.
4. Check4Cases(): This function will go through the 4 possible cases (JohnsonSU, shifted lognormal, shifted negative lognormal, constant) and find which is the best-fitting distribution to the given dataset, via least squares optimisation.
}
\value{
The following items will be returned no matter the optimal distribution:
\item{CaseID}{The distribution case in which the given dataset falls; one of JSU, Lognormal, Negative Lognormal or Constant.}
\item{SSQ}{The SSQ associated with the optimal distribution for the given dataset.}
\item{flags}{The values of flag1 and flag2. If flag1 is raised during the algorithm, it signifies multiple local minima may be present within the GSS. If flag2 is raised during the algorithm, it signifies that the fitted distribution has zero density at some observation(s) in the given dataset. Whether these flags indicate problems or not is at the discretion of the user.}

If the optimal distribution is the JSU distribution, the shifted lognormal distribution or the shifted negative lognormal distribution, the following will be returned aswell as CaseID, SSQ and flags:
\item{delta}{The value of delta associated with the optimal distribution.}
\item{xi}{The value of xi associated with the optimal distribution.}
\item{gamma}{The value of gamma associated with the optimal distribution.}
\item{lambda}{The value of lambda associated with the optimal distribution.}

If the optimal distribution is the constant distribution, the following will be returned aswell as CaseID, SSQ and flags:
\item{constant}{The value of constant associated with the optimal distribution. Note: this only occurs in the constant distribution.}

If the optimal distribution is the normal distribution, the following will be returned aswell as CaseID, SSQ and flags:
\item{mean}{The mean of the normal distribution.}
\item{std}{The standard deviation of the normal distribution.}

If the optimal distribution is the degenerate case, the following will be returned aswell as CaseID, SSQ and flags:
\item{smallest}{The smallest observation in the given dataset.}
\item{largest}{The largest observation in the given dataset.}
\item{average}{The average of the interior observations in the given dataset, i.e. excluding the smallest and largest.}
}
\references{
E. Cai,
Scripts and Functions: Using R to Implement the
Golden Section Search Method for Numerical Optimization, The
Chemical Statistician Blog, (2013), available at
https://chemicalstatistician.wordpress.com/2013/04/22/using-r-to-implement-the-golden-bisection-method/
}
\author{
%%  ~~who you are~~
CJ Clarke
}
\note{
%%  ~~further notes~~
Typically, users of this package need only use the JSUparameters() function directly and all other functions are used in the background.
}
\examples{
### simple example
dat = rnorm(500)
JSUparameters(dat)

### stock example (taken from my Master's thesis)
\dontrun{
# go to https://finance.yahoo.com/quote/%5EISEQ/history?period1=868838400&period2=1622419200&interval=1wk&filter=history&frequency=1wk&includeAdjustedClose=true
# download the csv file
# load in the csv file (check your own file path)
iseq = read.csv("Downloads/^ISEQ.csv")
# only store the closing price (convert from characters)
iseq = as.numeric(iseq$Close)
# interpolate the missing observation
which(is.na(iseq)) # 195
iseq[195] = (iseq[194] + iseq[196])/2
# calculate the log returns
returns = log(iseq[2:1247] / iseq[1:1246])
# find the best-fitting JSU distribution to this data
JSUparameters(returns)
}

### use your own stock returns data to find the best-fitting JSU distribution
}
