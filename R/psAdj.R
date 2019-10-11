#' @export
#' @importFrom survey svydesign
#' @importFrom tableone CreateTableOne ExtractSmd
#' @importFrom Matching Match

psAdj <- function(data,methods='all'){

  adjResults <- list()
  # genData <- data
  library("survey")
  library("tableone")
  library ("Matching")
  cNames <- colnames(genData)[-c(1:4)]
  pFormula <- paste(cNames,collapse=" + ")
  pFormula <- as.formula(paste("exposure"," ~ ",pFormula,sep=""))
  ps <- glm(pFormula, data=genData, family="binomial")
  genData$ps  = fitted(ps) ## V1-v25 were used for PS model ##
  genData$w.att<-ifelse(genData$exposure==1,1,genData$ps/(1-genData$ps));range(genData$w.att)  ## ATT weights in full data ##

  adjResults$propensity.score <- genData$ps

  ##### PS matching #####
  rr2 <- Match(Y=genData$obs.outcome,Tr=genData$exposure,X=genData$ps,M=1,replace=F,caliper=0.2); ## Match on PS score with 1:1 ratio
  matched.treated = unique(genData[rr2$index.treated,])
  matched.treated$Id = 1:length(matched.treated$ps)
  matched.untreated = genData[rr2$index.control,]
  matched.untreated$Id = rep(1:length(matched.treated$ps),each = 1) ## Ratio is 1:1
  matched_genDataa <- rbind(matched.treated,matched.untreated) ## Matched genDataaset


  ### precedent step for svyglm: combines a genDataa frame and all the survey design information ###
  design.att <- svydesign(ids=~1,weights= ~w.att,data=genData)
  #design.att2 <- svydesign(ids=~1,weights= ~w.att,data=genDataatt)

  vars <- colnames(genData)[5:ncol(genData)]
  tabraw <- CreateTableOne(vars = vars, strata = "exposure", data = genData, test = TRUE)
  smd.unadj = ExtractSmd(tabraw)

  tabweighted <- svyCreateTableOne(vars = vars, strata = "exposure", data = design.att, test = TRUE);
  smd.weighted = ExtractSmd(tabweighted)

  tabMatched <- CreateTableOne(vars = vars, strata = "exposure", data = matched_genDataa, test = TRUE)
  smd.matched = ExtractSmd(tabMatched)

  ## main analysis ###
  model0 <- glm(obs.outcome ~ exposure, family = "binomial", data=genData) ## Unadjusted
  model1 <- glm(obs.outcome ~ exposure, family = "quasibinomial", weights = w.att, data=genData)
  #model2 <- glm(obs.outcome ~ exposure + x1 + x2+ x3+ x4+x5+x6+x7+ x8+ x9+ x10 +x11+x12+x13+x14+x15+x16+x17+x18+x19+x20, family = "binomial", data=genData)
  model3 <- glm(obs.outcome ~ exposure, family = "binomial", data=matched_genDataa)

  # summary(model3)
  # round((summary(model3)$coef)[2,],5)

  naive = (summary(model0)$coef)[2,]
  IPW = (summary(model1)$coef)[2,]
  #est3 = (summary(model2)$coef)[2,]
  Matched = (summary(model3)$coef)[2,]

  smd.all = cbind(smd.unadj,smd.weighted,smd.matched,0)
  estimates <- rbind(naive,IPW,Matched)

  adjResults$estimates <- estimates

  unadj <-  data.frame(SMD = as.numeric(smd.unadj[(1:(dim(smd.unadj)[1]-2))]), Variable = rownames(smd.unadj)[1:(dim(smd.unadj)[1]-2)], Method = 'Unadjusted')
  weighted <- data.frame(SMD = as.numeric(smd.weighted[1:(dim(smd.weighted)[1]-2)]), Variable = rownames(smd.weighted)[1:(dim(smd.weighted)[1]-2)], Method = 'Weighted')
  matched <- data.frame(SMD = as.numeric(smd.matched[1:(dim(smd.matched)[1]-2)]), Variable = rownames(smd.matched)[1:(dim(smd.matched)[1]-2)], Method = 'Matched')
  smds <- rbind(unadj,weighted,matched)

  adjResults$SMD <- smds
  adjResults$data <- genData
  adjResults$matched.data <- matched_genDataa
  adjResults
}
