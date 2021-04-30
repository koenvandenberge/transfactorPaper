glmGamPoiResultsAll <- function(counts, design){
  library(glmGamPoi)
  # test all cell types
  res <- glm_gp(counts, design = design)
  ict <- rep(1, nrow(design))
  reduced_design <- model.matrix(~ -1 + ict)
  inf <- test_de(res, reduced_design = reduced_design)
  # for(kk in 1:(ncol(design)-1)){
  #   inf <- test_de(res, contrast = )
  # }
  #
  return(inf)
}

edgeRResultsAll <- function(counts, ct){
  library(edgeR)
  edgeDesign <- model.matrix(~ct)
  d <- DGEList(counts)
  d <- calcNormFactors(d)
  d <- estimateDisp(d, edgeDesign)
  fit <- glmFit(d, edgeDesign)
  lrt <- glmLRT(fit, coef=2:ncol(edgeDesign))
  return(lrt)
}

glmGamPoiResultsInd <- function(counts, design){
  library(glmGamPoi)
  # test each cell type individually vs reference
  res <- glm_gp(counts, design = design)
  contrast <- rep(0, ncol(design))
  resList <- list()
  for(kk in 1:(ncol(design)-1)){
    curContrast <- contrast
    curContrast[c(1,kk+1)] <- c(-1,1)
    inf <- test_de(res, contrast = curContrast)
    resList[[kk]] <- inf
  }
  return(resList)
}

limmaResultsInd <- function(data, design){
  require(limma)
  # test each cell type individually vs reference
  fits <- limma::lmFit(data, design)
  contrast <- rep(0, ncol(design))
  resList <- list()
  for(kk in 1:(ncol(design)-1)){
    curContrast <- contrast
    curContrast[c(1,kk+1)] <- c(-1,1)
    contFit <- contrasts.fit(fits, contrasts = curContrast)
    contFit <- eBayes(contFit)
    tt <- topTable(contFit, n=Inf, sort.by="none")
    resList[[kk]] <- tt
  }
  return(resList)
}

# limmaResultsAll <- function(data, design){
#   require(limma)
#   # test each cell type individually vs reference
#   fits <- limma::lmFit(data, design)
#   contrast.matrix <- makeContrasts(ctb-cta, ctc-cta, ctd-cta, cte-cta,
#                                    levels=design)
#   contFit <- contrasts.fit(fits, contrasts = contrast.matrix)
#   contFit <- eBayes(contFit)
#   tt <- topTableF(contFit, n=Inf, sort.by="none")
#   return(tt)
# }


limmaResultsAll <- function(data, design){
  require(limma)
  # test each cell type individually vs reference
  fits <- limma::lmFit(data, design)
  contrast.matrix <- matrix(0, nrow=ncol(design), ncol=ncol(design)-1,
                            dimnames = list(colnames(design),
                                            paste0(colnames(design)[-1],"-",colnames(design)[1])))
  for(cc in 1:ncol(contrast.matrix)) contrast.matrix[c(cc+1, 1),cc] <- c(1,-1)
  contFit <- contrasts.fit(fits, contrasts = contrast.matrix)
  contFit <- eBayes(contFit)
  tt <- topTableF(contFit, n=Inf, sort.by="none")
  return(tt)
}


simple_auc <- function(TPR, FPR){
  # function from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html
  # inputs already sorted, best scores first
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}


addNoiseToGRN <- function(X,
                          nNoiseTFs = 0,
                          nNoiseEdges = 1000,
                          noiseProb = 3/nrow(X),
                          seed = NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }

  XNoisy <- X

  ## add false edges
  if(nNoiseEdges > 0){
    # add 1000 noisy edges
    noiseEdges <- sample(which(X==0), size=nNoiseEdges)
    XNoisy[noiseEdges] <- 1
  } else if(nNoiseEdges < 0){
    allTrueEdges <- which(X==1)
    rmTrueEdges <- sample(allTrueEdges, abs(nNoiseEdges))
    XNoisy[rmTrueEdges] <- 0
  }

  ## add false TFs
  if(nNoiseTFs > 0){
    noiseTFX <- matrix(rbinom(n=nrow(X) * nNoiseTFs, size=1, p=noiseProb),
                       nrow=nrow(X), ncol=nNoiseTFs)
    while(!all(colSums(noiseTFX) > 0)){ #guarantee at least 1 gene
      id0 <- which(colSums(noiseTFX) == 0)
      noiseTFX[,id0] <- matrix(rbinom(n=nrow(X) * length(id0), size=1, p=noiseProb),
                               nrow=nrow(X), ncol=length(id0))
    }
    colnames(noiseTFX) <- paste0("noiseTF",1:ncol(noiseTFX))
    XNoisy <- cbind(XNoisy, noiseTFX)
  }

  return(XNoisy)
}

constructViperRegulon <- function(X, alpha){
  tfAll <- unlist(mapply(rep, colnames(X), each=colSums(abs(X))))
  targetAll <- rownames(X)[unlist(apply(abs(X),2, function(x) which(x > 0)))]
  mor <- X[abs(X)>0]
  alphaAll <- alpha[abs(X)>0] / max(alpha)
  dfReg <- data.frame(tf=tfAll,
                      target=targetAll,
                      mor=X[abs(X)>0],
                      likelihood=alphaAll)
  dfReg <- dfReg[!duplicated(dfReg),]
  dfReg$tf <- as.character(dfReg$tf)
  dfReg$target <- as.character(dfReg$target)
  regulon <- dorothea:::df2regulon(dfReg)
  return(regulon)
}

constructGenesets <- function(X, alpha){
  tfAll <- unlist(mapply(rep, colnames(X), each=colSums(abs(X))))
  targetAll <- rownames(X)[unlist(apply(abs(X),2, function(x) which(x > 0)))]
  mor <- X[abs(X)>0]
  alphaAll <- alpha[abs(X)>0] / max(alpha)
  dfReg <- data.frame(tf=tfAll,
                      target=targetAll,
                      mor=X[abs(X)>0],
                      likelihood=alphaAll)
  dfReg <- dfReg[!duplicated(dfReg),]
  dfReg$tf <- as.character(dfReg$tf)
  dfReg$target <- as.character(dfReg$target)
  genesets = dfReg %>%
    group_by(tf) %>%
    summarise(geneset = list(GSEABase::GeneSet(target))) %>%
    transmute(tf, geneset2 = pmap(., .f=function(tf, geneset, ...) {
      setName(geneset) = tf
      return(geneset)
    })) %>%
    deframe() %>%
    GeneSetCollection()
}

evaluateSimulation_noRepressions <- function(counts,
                               design,
                               X,
                               alpha,
                               regulon,
                               genesets,
                               truth,
                               verbose = FALSE,
                               alphaScale = 1){


  ## Poisson model
  poisRes <- transfactor::estimateActivity(counts = counts,
                                      X = X,
                                      model = "poisson",
                                      U = design,
                                      verbose = verbose,
                                      plot = FALSE,
                                      maxIter = 500,
                                      epsilon = .1,
                                      sparse = FALSE,
                                      repressions = FALSE)
  Y_ti_pois <- transfactor::tfCounts(poisRes$mu_gtc, counts, design)
  fTestGamPoi_pois <- glmGamPoiResultsAll(Y_ti_pois, design)
  indTestGamPoi_pois <- glmGamPoiResultsInd(Y_ti_pois, design)

  ## Dir-Mult model, known alpha
  dirMultRes_alphaKnown <- transfactor::estimateActivity(counts = counts,
                                                 X = X,
                                                 model = "dirMult",
                                                 alpha = alpha,
                                                 U = design,
                                                 verbose = verbose,
                                                 plot = FALSE,
                                                 maxIter = 500,
                                                 epsilon=.1,
                                                 alphaScale = alphaScale,
                                                 sparse = FALSE,
                                                 repressions = FALSE)
  Y_ti_dirMult <- transfactor::tfCounts(dirMultRes_alphaKnown$mu_gtc, counts, design)
  fTestGamPoi_dirMult <- glmGamPoiResultsAll(Y_ti_dirMult, design)
  indTestGamPoi_dirMult <- glmGamPoiResultsInd(Y_ti_dirMult, design)

  ## Dir-Mult model, estimate alpha
  dirMultRes_alphaEst <- transfactor::estimateActivity(counts = counts,
                                                       X = X,
                                                       model = "dirMultAlpha",
                                                       U = design,
                                                       verbose = verbose,
                                                       plot = FALSE,
                                                       maxIter = 500,
                                                       epsilon=.1,
                                                       alphaScale = alphaScale,
                                                       sparse = FALSE,
                                                       repressions = FALSE)
  Y_ti_dirMult_alphaEst <- transfactor::tfCounts(dirMultRes_alphaEst$mu_gtc, counts, design)
  fTestGamPoi_dirMult_alphaEst <- glmGamPoiResultsAll(Y_ti_dirMult_alphaEst, design)
  indTestGamPoi_dirMult_alphaEst <- glmGamPoiResultsInd(Y_ti_dirMult_alphaEst, design)

  ## Poisson model, lasso
  poisResLasso <- transfactor::estimateActivity(counts = counts,
                                                X = X,
                                                model = "poisson",
                                                U = design,
                                                verbose = verbose,
                                                plot = FALSE,
                                                maxIter = 500,
                                                epsilon=.1,
                                                sparse = TRUE,
                                                repressions = FALSE)
  Y_ti_poisLasso <- transfactor::tfCounts(poisResLasso$mu_gtc, counts, design)
  fTestGamPoi_poisLasso <- glmGamPoiResultsAll(Y_ti_poisLasso, design)
  indTestGamPoi_poisLasso <- glmGamPoiResultsInd(Y_ti_poisLasso, design)

  ## Dir-Mult model, known alpha, lasso
  dirMultRes_alphaKnown_lasso <- transfactor::estimateActivity(counts = counts,
                                                               X = X,
                                                               model = "dirMult",
                                                               alpha = alpha,
                                                               U = design,
                                                               verbose = verbose,
                                                               plot = FALSE,
                                                               maxIter = 500,
                                                               epsilon=.1,
                                                               alphaScale = alphaScale,
                                                               sparse = TRUE,
                                                               repressions = FALSE)
  Y_ti_dirMultLasso <- transfactor::tfCounts(dirMultRes_alphaKnown_lasso$mu_gtc, counts, design)
  fTestGamPoi_dirMultLasso <- glmGamPoiResultsAll(Y_ti_dirMultLasso, design)
  indTestGamPoi_dirMultLasso <- glmGamPoiResultsInd(Y_ti_dirMultLasso, design)

  ## Dir-Mult model, estimate alpha, lasso
  dirMultRes_alphaEst_lasso <- transfactor::estimateActivity(counts = counts,
                                                              X = X,
                                                             model = "dirMultAlpha",
                                                              U = design,
                                                              verbose = verbose,
                                                              plot = FALSE,
                                                              maxIter = 500,
                                                              epsilon = .1,
                                                              sparse = TRUE,
                                                              repressions = FALSE)
  Y_ti_dirMult_alphaEst_lasso <- transfactor::tfCounts(dirMultRes_alphaEst_lasso$mu_gtc, counts, design)
  fTestGamPoi_dirMult_alphaEst_lasso <- glmGamPoiResultsAll(Y_ti_dirMult_alphaEst_lasso, design)
  indTestGamPoi_dirMult_alphaEst_lasso <- glmGamPoiResultsInd(Y_ti_dirMult_alphaEst_lasso, design)

  ## viper
  viperRes <- viper(counts,
                    regulon,
                    nes = TRUE,
                    method = "scale",
                    minsize = 2,
                    eset.filter = F,
                    adaptive.size = F)
  viperRes <- viperRes[colnames(X)[colnames(X) %in% rownames(viperRes)],]
  indTestLimma_viper <- limmaResultsInd(viperRes, design)
  fTestLimma_viper <- limmaResultsAll(viperRes, design)

  ## AUCell
  obj <- AUCell_buildRankings(t(scale(t(counts))), nCores=1, plotStats = F, verbose = F) %>%
    AUCell_calcAUC(genesets, ., verbose=F)
  resAUCell <- AUCell::getAUC(obj)
  resAUCell <- resAUCell[colnames(X),]
  indTestLimma_AUCell <- limmaResultsInd(resAUCell, design)
  fTestLimma_AUCell <- limmaResultsAll(resAUCell, design)

  ## ROC curves for each cell type
  plistInd <- list()
  for(kk in 1:(ncol(design)-1)){
    # DE id
    deId <- truth[,kk]

    # poisson
    ooF <- order(indTestGamPoi_pois[[kk]]$f_statistic, decreasing = TRUE)
    tpr <- cumsum(deId[ooF]) / sum(deId)
    fpr <- cumsum(!deId[ooF]) / sum(!deId)
    df <- data.frame(tpr=tpr,
                     fpr=fpr,
                     method="poisson")

    # poisson, gaussian lasso
    ooF <- order(indTestGamPoi_poisLasso[[kk]]$f_statistic, decreasing = TRUE)
    tpr <- cumsum(deId[ooF]) / sum(deId)
    fpr <- cumsum(!deId[ooF]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="poisson_lasso")
    df <- rbind(df, curDf)


    # dirMult
    ooF <- order(indTestGamPoi_dirMult[[kk]]$f_statistic, decreasing = TRUE)
    tpr <- cumsum(deId[ooF]) / sum(deId)
    fpr <- cumsum(!deId[ooF]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="dirMult")
    df <- rbind(df, curDf)

    # dirMult, gaussian lasso
    ooF <- order(indTestGamPoi_dirMultLasso[[kk]]$f_statistic, decreasing = TRUE)
    tpr <- cumsum(deId[ooF]) / sum(deId)
    fpr <- cumsum(!deId[ooF]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="dirMult_lasso")
    df <- rbind(df, curDf)

    # dirMult, est alpha
    ooF <- order(indTestGamPoi_dirMult_alphaEst[[kk]]$f_statistic, decreasing = TRUE)
    tpr <- cumsum(deId[ooF]) / sum(deId)
    fpr <- cumsum(!deId[ooF]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="dirMultEstAlpha")
    df <- rbind(df, curDf)

    # dirMult, est alpha, lasso
    ooF <- order(indTestGamPoi_dirMult_alphaEst_lasso[[kk]]$f_statistic, decreasing = TRUE)
    tpr <- cumsum(deId[ooF]) / sum(deId)
    fpr <- cumsum(!deId[ooF]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="dirMultEstAlpha_lasso")
    df <- rbind(df, curDf)


    # viper
    ooVip <- order(abs(indTestLimma_viper[[kk]]$t), decreasing = TRUE)
    tpr <- cumsum(deId[ooVip]) / sum(deId)
    fpr <- cumsum(!deId[ooVip]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="viper")
    df <- rbind(df, curDf)

    # AUCell
    ooAUC <- order(abs(indTestLimma_AUCell[[kk]]$t), decreasing = TRUE)
    tpr <- cumsum(deId[ooAUC]) / sum(deId)
    fpr <- cumsum(!deId[ooAUC]) / sum(!deId)
    curDf <- data.frame(tpr=tpr,
                        fpr=fpr,
                        method="AUCell")
    df <- rbind(df, curDf)

    df$ct <- kk

    if(kk == 1){
      dfIndAll <- df
    } else {
      dfIndAll <- rbind(dfIndAll, df)
    }

    plistInd[[kk]] <- ggplot(df, aes(x=fpr, y=tpr, gorup=method, col=method)) +
      geom_path() +
      theme_classic()
  }


  ## ROC curves across all cell types
  deId <- truth[,1]

  # poisson
  ooF <- order(fTestGamPoi_pois$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  df <- data.frame(tpr=tpr,
                   fpr=fpr,
                   method="poisson")
  auc_pois <- simple_auc(tpr, fpr)

  # poisson,  lasso
  ooF <- order(fTestGamPoi_poisLasso$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="poisson_lasso")
  df <- rbind(df, curDf)
  auc_poisLasso <- simple_auc(tpr, fpr)


  # dirMult
  ooF <- order(fTestGamPoi_dirMult$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMult")
  df <- rbind(df, curDf)
  auc_dirMult <- simple_auc(tpr, fpr)

  # dirMult,  lasso
  ooF <- order(fTestGamPoi_dirMultLasso$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMult_lasso")
  df <- rbind(df, curDf)
  auc_dirMultLasso <- simple_auc(tpr, fpr)

  # dirMult, est alpha
  ooF <- order(fTestGamPoi_dirMult_alphaEst$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMultEstAlpha")
  df <- rbind(df, curDf)
  auc_dirMult_alphaEst <- simple_auc(tpr, fpr)

  # dirMult, est alpha, lasso
  ooF <- order(fTestGamPoi_dirMult_alphaEst_lasso$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMultEstAlpha_lasso")
  df <- rbind(df, curDf)
  auc_dirMult_alphaEst_lasso <- simple_auc(tpr, fpr)

  # viper
  ooVip <- order(abs(fTestLimma_viper$F), decreasing = TRUE)
  tpr <- cumsum(deId[ooVip]) / sum(deId)
  fpr <- cumsum(!deId[ooVip]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="viper")
  df <- rbind(df, curDf)
  auc_viper <- simple_auc(tpr, fpr)

  # AUCell
  ooAUC <- order(abs(fTestLimma_AUCell$F), decreasing = TRUE)
  tpr <- cumsum(deId[ooAUC]) / sum(deId)
  fpr <- cumsum(!deId[ooAUC]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="AUCell")
  df <- rbind(df, curDf)
  auc_AUCell <- simple_auc(tpr, fpr)

  pAll <- ggplot(df, aes(x=fpr, y=tpr, gorup=method, col=method)) +
    geom_path() +
    theme_classic()

  dfAUC <- data.frame(auc = c(auc_pois,
                              auc_dirMult,
                              auc_dirMult_alphaEst,
                              auc_poisLasso,
                              auc_dirMultLasso,
                              auc_dirMult_alphaEst_lasso,
                              auc_viper,
                              auc_AUCell),
                      method = c("poisson",
                                 "dirMult",
                                 "dirMult_alphaEst",
                                 "poisson_lasso",
                                 "dirMult_lasso",
                                 "dirMult_lasso_alphaEst",
                                 "viper",
                                 "AUCell"))


  return(list(plistInd = plistInd,
              dfIndAll = dfIndAll,
              dfAll = df,
              pAll = pAll,
              dfAUC = dfAUC))
}


# counts = countsAll
# design = design
# X = XNoisy
# alpha = alphaTrueNoisy
# regulon = regulonNoisy
# genesets = genesetsNoisy
# truth = truthNoisy
# verbose = TRUE
# alphaScale = 1
# iterMax = 500

evaluateSimulation_repressions <- function(counts,
                                            design,
                                            X,
                                            alpha,
                                            regulon,
                                            genesets,
                                            truth,
                                            verbose = FALSE,
                                            alphaScale = 1,
                                           iterMax = 500){

  XPos <- X
  XPos[XPos == -1] <- 0
  alphaPos <- alpha
  alphaPos[XPos == 0] <- 0
  keepXPosRow <- rowSums(XPos) > 0
  XPos <- XPos[keepXPosRow,]
  alphaPos <- alphaPos[keepXPosRow,]
  keepXPosCol <- colSums(XPos)>0
  XPos <- XPos[,keepXPosCol]
  alphaPos <- alphaPos[,keepXPosCol]
  countsPos <- counts[rownames(XPos),]


  ## Poisson model
  poisRes <- emPoissonThinning_faster(counts = countsPos,
                                      X = XPos,
                                      rho_t = NULL,
                                      U = design,
                                      verbose = verbose,
                                      plot = FALSE,
                                      nIters = iterMax,
                                      epsilon=.1)
  Y_ti_pois <- getCellLevelYt_faster(poisRes$mu_gti, countsPos)
  fTestGamPoi_pois <- glmGamPoiResultsAll(Y_ti_pois, design)
  indTestGamPoi_pois <- glmGamPoiResultsInd(Y_ti_pois, design)

  ## Dir-Mult model, known alpha
  dirMultRes_alphaKnown <- emDirMult_alphaKnown3(counts = countsPos,
                                                 X = XPos,
                                                 alpha = alphaPos,
                                                 rho_t = NULL,
                                                 U = design,
                                                 verbose = verbose,
                                                 plot = FALSE,
                                                 nIters = iterMax,
                                                 epsilon=.1,
                                                 alphaScale = alphaScale)
  Y_ti_dirMult <- getCellLevelYt_faster(dirMultRes_alphaKnown$mu_gti, countsPos)
  fTestGamPoi_dirMult <- glmGamPoiResultsAll(Y_ti_dirMult, design)
  indTestGamPoi_dirMult <- glmGamPoiResultsInd(Y_ti_dirMult, design)

  ## Dir-Mult model, estimate alpha
  dirMultRes_alphaEst <- emDirMult_estAlpha_iterative2(counts = countsPos,
                                                       X = XPos,
                                                       rho_t = NULL,
                                                       U = design,
                                                       verbose = verbose,
                                                       plot = FALSE,
                                                       nIters = iterMax,
                                                       epsilon=.1)
  Y_ti_dirMult_alphaEst <- getCellLevelYt_faster(dirMultRes_alphaEst$mu_gti, countsPos)
  fTestGamPoi_dirMult_alphaEst <- glmGamPoiResultsAll(Y_ti_dirMult_alphaEst, design)
  indTestGamPoi_dirMult_alphaEst <- glmGamPoiResultsInd(Y_ti_dirMult_alphaEst, design)

  ## Poisson model, lasso
  poisResLasso <- emPoissonThinning_faster_lassoInit_sufStat(counts = countsPos,
                                                             X = XPos,
                                                             rho_t = NULL,
                                                             U = design,
                                                             verbose = verbose,
                                                             plot = FALSE,
                                                             nIters = iterMax,
                                                             epsilon=.1)
  Y_ti_poisLasso <- getCellLevelYti_faster_sufStats(poisResLasso$mu_gtc, countsPos, design)
  fTestGamPoi_poisLasso <- glmGamPoiResultsAll(Y_ti_poisLasso, design)
  indTestGamPoi_poisLasso <- glmGamPoiResultsInd(Y_ti_poisLasso, design)

  ## Dir-Mult model, known alpha, lasso
  dirMultRes_alphaKnown_lasso <- emDirMultThinning_faster_lassoInit_sufStat(counts = countsPos,
                                                                            X = XPos,
                                                                            alpha = alphaPos,
                                                                            rho_t = NULL,
                                                                            U = design,
                                                                            verbose = verbose,
                                                                            plot = FALSE,
                                                                            nIters = iterMax,
                                                                            epsilon=.1,
                                                                            alphaScale = alphaScale)
  Y_ti_dirMult_lasso <- getCellLevelYti_faster_sufStats(dirMultRes_alphaKnown_lasso$mu_gtc, countsPos, design)
  fTestGamPoi_dirMultLasso <- glmGamPoiResultsAll(Y_ti_dirMult_lasso, design)
  indTestGamPoi_dirMultLasso <- glmGamPoiResultsInd(Y_ti_dirMult_lasso, design)

  ## Dir-Mult model, estimate alpha, lasso
  dirMultRes_alphaEst_lasso <- emDirMultThinning_faster_lassoInit_sufStat_repressions_estAlpha3(
    counts = countsPos,
    X = XPos,
    rho_t = NULL,
    U = design,
    verbose = verbose,
    plot = FALSE,
    nIters = iterMax,
    epsilon = .1,
    repressions = FALSE)
  Y_ti_dirMult_alphaEst_lasso <- getCellLevelYti_faster_sufStats(dirMultRes_alphaEst_lasso$mu_gtc, countsPos, design)
  fTestGamPoi_dirMult_alphaEst_lasso <- glmGamPoiResultsAll(Y_ti_dirMult_alphaEst_lasso, design)
  indTestGamPoi_dirMult_alphaEst_lasso <- glmGamPoiResultsInd(Y_ti_dirMult_alphaEst_lasso, design)

  ## Poisson model, lasso, repressions
  poisResLasso_repr <- emPoissonThinning_faster_lassoInit_sufStat_repressions(counts = counts,
                                                             X = X,
                                                             rho_t = NULL,
                                                             U = design,
                                                             verbose = verbose,
                                                             plot = FALSE,
                                                             nIters = iterMax,
                                                             epsilon=.1,
                                                             repressions = TRUE)
  Y_ti_poisLasso_repr <- getCellLevelYti_faster_sufStats(poisResLasso_repr$mu_gtc, counts, design)
  fTestGamPoi_poisLasso_repr <- glmGamPoiResultsAll(Y_ti_poisLasso_repr, design)
  indTestGamPoi_poisLasso_repr <- glmGamPoiResultsInd(Y_ti_poisLasso_repr, design)

  ## Dir-Mult model, known alpha, lasso, repressions
  dirMultRes_alphaKnown_lasso_repr <- emDirMultThinning_faster_lassoInit_sufStat_repressions(counts = counts,
                                                                            X = X,
                                                                            alpha = alpha,
                                                                            rho_t = NULL,
                                                                            U = design,
                                                                            verbose = verbose,
                                                                            plot = FALSE,
                                                                            nIters = iterMax,
                                                                            epsilon=.1,
                                                                            alphaScale = alphaScale,
                                                                            repressions = TRUE)
  Y_ti_dirMultLasso_repr <- getCellLevelYti_faster_sufStats(dirMultRes_alphaKnown_lasso_repr$mu_gtc, counts, design)
  fTestGamPoi_dirMultLasso_repr <- glmGamPoiResultsAll(Y_ti_dirMultLasso_repr, design)
  indTestGamPoi_dirMultLasso_repr <- glmGamPoiResultsInd(Y_ti_dirMultLasso_repr, design)

  ## Dir-Mult model, estimate alpha, lasso, repressions
  dirMultRes_alphaEst_lasso_repr <- emDirMultThinning_faster_lassoInit_sufStat_repressions_estAlpha3(
    counts = counts,
    X = X,
    rho_t = NULL,
    U = design,
    verbose = verbose,
    plot = FALSE,
    nIters = iterMax,
    epsilon = .1,
    repressions = TRUE)
  Y_ti_dirMult_alphaEst_lasso_repr <- getCellLevelYti_faster_sufStats(dirMultRes_alphaEst_lasso_repr$mu_gtc, counts, design)
  fTestGamPoi_dirMult_alphaEst_lasso_repr <- glmGamPoiResultsAll(Y_ti_dirMult_alphaEst_lasso_repr, design)
  indTestGamPoi_dirMult_alphaEst_lasso_repr <- glmGamPoiResultsInd(Y_ti_dirMult_alphaEst_lasso_repr, design)

  ## viper
  viperRes <- viper(counts,
                    regulon,
                    nes = TRUE,
                    method = "scale",
                    minsize = 2,
                    eset.filter = F,
                    adaptive.size = F)
  viperRes <- viperRes[colnames(X)[colnames(X) %in% rownames(viperRes)],]
  indTestLimma_viper <- limmaResultsInd(viperRes, design)
  fTestLimma_viper <- limmaResultsAll(viperRes, design)

  ## AUCell
  obj <- AUCell_buildRankings(t(scale(t(counts))), nCores=1, plotStats = F, verbose = F) %>%
    AUCell_calcAUC(genesets, ., verbose=F)
  resAUCell <- AUCell::getAUC(obj)
  resAUCell <- resAUCell[colnames(X),]
  indTestLimma_AUCell <- limmaResultsInd(resAUCell, design)
  fTestLimma_AUCell <- limmaResultsAll(resAUCell, design)

  # ## ROC curves for each cell type
  # plistInd <- list()
  # for(kk in 1:(ncol(design)-1)){
  #   # DE id
  #   deId <- truth[,kk]
  #
  #   # poisson
  #   ooF <- order(indTestGamPoi_pois[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   df <- data.frame(tpr=tpr,
  #                    fpr=fpr,
  #                    method="poisson")
  #
  #   # poisson, lasso
  #   ooF <- order(indTestGamPoi_poisLasso[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="poisson_lasso")
  #   df <- rbind(df, curDf)
  #
  #   # poisson, lasso, repressions
  #   ooF <- order(indTestGamPoi_poisLasso_repr[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="poisson_lasso_repr")
  #   df <- rbind(df, curDf)
  #
  #
  #   # dirMult
  #   ooF <- order(indTestGamPoi_dirMult[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="dirMult")
  #   df <- rbind(df, curDf)
  #
  #   # dirMult, lasso
  #   ooF <- order(indTestGamPoi_dirMultLasso[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="dirMult_lasso")
  #   df <- rbind(df, curDf)
  #
  #   # dirMult, lasso, repressions
  #   ooF <- order(indTestGamPoi_dirMultLasso_repr[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="dirMult_lasso_repr")
  #   df <- rbind(df, curDf)
  #
  #   # dirMult, est alpha
  #   ooF <- order(indTestGamPoi_dirMult_alphaEst[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="dirMultEstAlpha")
  #   df <- rbind(df, curDf)
  #
  #   # dirMult, est alpha, lasso
  #   ooF <- order(indTestGamPoi_dirMult_alphaEst_lasso[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="dirMultEstAlpha_lasso")
  #   df <- rbind(df, curDf)
  #
  #   # dirMult, est alpha, lasso, repressions
  #   ooF <- order(indTestGamPoi_dirMult_alphaEst_lasso_repr[[kk]]$f_statistic, decreasing = TRUE)
  #   tpr <- cumsum(deId[ooF]) / sum(deId)
  #   fpr <- cumsum(!deId[ooF]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="dirMultEstAlpha_lasso_repr")
  #   df <- rbind(df, curDf)
  #
  #
  #   # viper
  #   ooVip <- order(abs(indTestLimma_viper[[kk]]$t), decreasing = TRUE)
  #   tpr <- cumsum(deId[ooVip]) / sum(deId)
  #   fpr <- cumsum(!deId[ooVip]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="viper")
  #   df <- rbind(df, curDf)
  #
  #   # AUCell
  #   ooAUC <- order(abs(indTestLimma_AUCell[[kk]]$t), decreasing = TRUE)
  #   tpr <- cumsum(deId[ooAUC]) / sum(deId)
  #   fpr <- cumsum(!deId[ooAUC]) / sum(!deId)
  #   curDf <- data.frame(tpr=tpr,
  #                       fpr=fpr,
  #                       method="AUCell")
  #   df <- rbind(df, curDf)
  #
  #   df$ct <- kk
  #
  #   if(kk == 1){
  #     dfIndAll <- df
  #   } else {
  #     dfIndAll <- rbind(dfIndAll, df)
  #   }
  #
  #   plistInd[[kk]] <- ggplot(df, aes(x=fpr, y=tpr, gorup=method, col=method)) +
  #     geom_path() +
  #     theme_classic()
  # }


  ## ROC curves across all cell types
  deId <- truth[,1]

  # poisson
  ooF <- order(fTestGamPoi_pois$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  df <- data.frame(tpr=tpr,
                   fpr=fpr,
                   method="poisson")
  auc_pois <- simple_auc(tpr, fpr)

  # poisson,  lasso
  ooF <- order(fTestGamPoi_poisLasso$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="poisson_lasso")
  df <- rbind(df, curDf)
  auc_poisLasso <- simple_auc(tpr, fpr)

  # poisson,  lasso, repressions
  ooF <- order(fTestGamPoi_poisLasso_repr$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="poisson_lasso_repr")
  df <- rbind(df, curDf)
  auc_poisLasso_repr <- simple_auc(tpr, fpr)

  # dirMult
  ooF <- order(fTestGamPoi_dirMult$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMult")
  df <- rbind(df, curDf)
  auc_dirMult <- simple_auc(tpr, fpr)

  # dirMult,  lasso
  ooF <- order(fTestGamPoi_dirMultLasso$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMult_lasso")
  df <- rbind(df, curDf)
  auc_dirMultLasso <- simple_auc(tpr, fpr)

  # dirMult,  lasso, repressions
  ooF <- order(fTestGamPoi_dirMultLasso_repr$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMult_lasso_repr")
  df <- rbind(df, curDf)
  auc_dirMultLasso_repr <- simple_auc(tpr, fpr)

  # dirMult, est alpha
  ooF <- order(fTestGamPoi_dirMult_alphaEst$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMultEstAlpha")
  df <- rbind(df, curDf)
  auc_dirMult_alphaEst <- simple_auc(tpr, fpr)

  # dirMult, est alpha, lasso
  ooF <- order(fTestGamPoi_dirMult_alphaEst_lasso$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMultEstAlpha_lasso")
  df <- rbind(df, curDf)
  auc_dirMult_alphaEst_lasso <- simple_auc(tpr, fpr)

  # dirMult, est alpha, lasso, repressions
  ooF <- order(fTestGamPoi_dirMult_alphaEst_lasso_repr$f_statistic, decreasing = TRUE)
  tpr <- cumsum(deId[ooF]) / sum(deId)
  fpr <- cumsum(!deId[ooF]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="dirMultEstAlpha_lasso_repr")
  df <- rbind(df, curDf)
  auc_dirMult_alphaEst_lasso_repr <- simple_auc(tpr, fpr)

  # viper
  ooVip <- order(abs(fTestLimma_viper$F), decreasing = TRUE)
  tpr <- cumsum(deId[ooVip]) / sum(deId)
  fpr <- cumsum(!deId[ooVip]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="viper")
  df <- rbind(df, curDf)
  auc_viper <- simple_auc(tpr, fpr)

  # AUCell
  ooAUC <- order(abs(fTestLimma_AUCell$F), decreasing = TRUE)
  tpr <- cumsum(deId[ooAUC]) / sum(deId)
  fpr <- cumsum(!deId[ooAUC]) / sum(!deId)
  curDf <- data.frame(tpr=tpr,
                      fpr=fpr,
                      method="AUCell")
  df <- rbind(df, curDf)
  auc_AUCell <- simple_auc(tpr, fpr)

  pAll <- ggplot(df, aes(x=fpr, y=tpr, gorup=method, col=method)) +
    geom_path() +
    theme_classic()

  dfAUC <- data.frame(auc = c(auc_pois,
                              auc_dirMult,
                              auc_dirMult_alphaEst,
                              auc_poisLasso,
                              auc_dirMultLasso,
                              auc_dirMult_alphaEst_lasso,
                              auc_poisLasso_repr,
                              auc_dirMultLasso_repr,
                              auc_dirMult_alphaEst_lasso_repr,
                              auc_viper,
                              auc_AUCell),
                      method = c("poisson",
                                 "dirMult",
                                 "dirMult_alphaEst",
                                 "poisson_lasso",
                                 "dirMult_lasso",
                                 "dirMult_lasso_alphaEst",
                                 "poisson_lasso_repr",
                                 "dirMult_lasso_repr",
                                 "dirMult_lasso_alphaEst_repr",
                                 "viper",
                                 "AUCell"))


  return(list(#plistInd = plistInd,
              #dfIndAll = dfIndAll,
              dfAll = df,
              pAll = pAll,
              dfAUC = dfAUC))
}

