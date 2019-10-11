\name{psAdj}
\alias{psAdj}

\title{
Propensity Score Adjustment
}
\description{
Performs propensity score adjustment
}
\usage{
psAdj(data,methods='all')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{}
    \item{methods}{}
}
\value{
 \item{SMD}{Standardized mean difference}
 \item{propensity.score}{Propensity score}
}
\author{
Jaejoon Song, Tae Hyun Jung
}
\examples{
# Load package
library(propensity)

# Simulating data with 20 confounders with sample size of 1000
genData <- simConf(n = 1000,nconf = 30,dconf = c('beta',2,2),
                   nexpred = 5, dexpred = c('beta',2,2),
                   noutpred = 5, doutpred = c('beta',2,2),
                   effect = 2, escale = 'OR',
                   brate = .1,erate =.1)

# Performing propensity score adjustment
adj <- psAdj(data=genData)

# Diagnostic plots
dPlot(adjResults=adj)
}
