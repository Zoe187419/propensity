\name{simConf}
\alias{simConf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulates data with confounders
}
\description{
Generates user specified confounded data to conduct simulation studies for causal inference methods
}
\usage{
simConf(n, nconf, dconf = c("beta", 2, 2),
           nexpred, dexpred = c("beta", 2, 2),
           noutpred, doutpred = c("beta", 2, 2),
           effect = log(2), escale = "OR", brate, erate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size}
  \item{nconf}{Number of confounders}
  \item{dconf}{Distribution of confounders}
  \item{nexpred}{Number of exposure predictors}
  \item{dexpred}{Distribution of exposure predictors}
  \item{noutpred}{Number of outcome predictors}
  \item{doutpred}{Distribution of outcome predictors}
  \item{effect}{Treatment effect}
  \item{escale}{Scale of treatment effect. Default is odds ratio.}
  \item{brate}{Outcome background rate}
  \item{erate}{Exposure rate}
}

\value{
A simulated dataset is returned
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
