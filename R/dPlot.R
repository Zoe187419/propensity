#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggpubr ggarrange
#' @importFrom plyr revalue


dPlot <- function(adjResults){

  #
  ## Plot
  library(ggplot2)

  # adjResults <- adj
  # adj$SMD

  fig = ggplot(data =  adjResults$SMD, mapping = aes(x = Variable, y = SMD, group = Method,color=Method, linetype = Method)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.1, size = 1,col="black",linetype=2) +
    coord_flip() +
    theme_bw() + theme(legend.key = element_blank()) + labs(y="Standardized Mean Difference")
  fig = fig +  facet_wrap(vars(Method))

  library(ggpubr)
  library(plyr)

  adjResults$data$exposure <- as.factor(adjResults$data$exposure)
  adjResults$data$exposure <- revalue(adjResults$data$exposure, c("1"="Exposed", "0"="Unexposed"))

  p1<-ggplot(adjResults$data, aes(x=ps, fill=as.factor(exposure))) +
    geom_density(alpha=0.4)+theme(legend.position="none") +
    scale_y_continuous(limits = c(0,6)) + ggtitle("Propensity Score (Unmatched Sample)")+ theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Propensity Score") + ylab("Density")

  p2<-ggplot(adjResults$matched.data, aes(x=ps, fill=as.factor(exposure))) +
    geom_density(alpha=0.4)+ scale_y_continuous(limits = c(0,6))  +
    ggtitle("Propensity Score (Unmatched Sample)") + theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Propensity Score") + ylab("Density")+ scale_fill_discrete(name = "Exposure")

  ggarrange(fig, p1, p2, nrow = 1)

}



