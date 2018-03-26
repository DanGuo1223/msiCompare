#' @title Fit hiearchical spatial model to MSI data
#' @description compareMSI is used to fit a hiearchical Bayesian spatial model to MSI data using a Gibbs Sampler MCMC approach. The model is fit separately for each m/z feature.
#' @param msset an object of class "MSImageSet"
#' @param conditionOfInterest a vector or factor giving the level of the condition of interest for each pixel in msset
#' @param feature the index of the m/z features for which the model should be fit
#' @param nsim number of desired MCMC samples
#' @param burnin number of MCMC samples to discard
#' @param trace logical, should the full list of MCMC samples be returned for each variable?
#' @param piPrior prior probability of differential abundance
#' @param seed random seed
#' @param logbase2 logical, should the intensities be log transformed?
#' @param coord data fram of coordinates of the MSImageSet, with columns 'x' and 'y'
#' @param type.neighbor neighborhood type (see adj.grid)
#' @param radius.neighbor desired neighborhood radius if neighborhood type 'radius' is selected (see adj.grid)
#' @param maxdist.neighbor maximum distance for locations to be considered neighbors if neighborhood type 'max.dist' is selected (see adj.grid)
#' @param spInit optional, provide precomputed spatial information from output of intializeSpatial
#' @param bioRep optional, vector or factor giving the individual/donor to which pixel in the msset belongs
#' @param techRep vector or factor giving the tissue to which each pixel in the msset belongs
#' @param beta0 prior mean of baseline effect
#' @param prec0 prior variance of baseline effect
#' @param precAlpha0 prior mean of condition 2 effect
#' @param a0_eps shape parameter for measurment error precision hyperprior
#' @param a0_bio shape parameter for biological replicate error precision hyperprior
#' @param b0_eps scale parameter for measurment error precision hyperprior
#' @param b0_bio scale parameter for biological replicate error precision hyperprior
#' @param a0_tec shape parameter for sample to sample error precision hyperprior
#' @param b0_tec scale parameter for sample to sample error precision hyperprior
#' @param a0_sp shape parameter for spatial precision hyperprior
#' @param b0_sp scale parameter for spatial precision hyperprior
#' @param rd ratio of spike variance to slab variance for condition 2 effect
#' @return res
#' @import mvtnorm
#' @import lme4
#' @import spam
#' @import coda
#' @export
#'

compareMSI <- function(msset,conditionOfInterest,
                       feature, nsim=5000, burnin = 2500, trace = T,
                       piPrior = .1, seed = 1, logbase2 = F, coord = NULL,
                       type.neighbor = "radius", radius.neighbor = 1, maxdist.neighbor = NULL,
                       spInit = NULL,
                       bioRep = NULL,
                       techRep,
                       beta0 = 0, # Prior Mean for beta, only allow intercept
                       prec0 = .01, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                       precAlpha0 = .01, #Prior Precision of slab (value of condition effect if it is not zero)
                       a0_eps=.001, b0_eps=.001,			# Hyperprior for tau (1/eps.var)
                       a0_bio=.001, b0_bio=.001,			# Hyperprior for taubio
                       a0_tec=.001, b0_tec=.001,			# Hyperprior for tautec
                       a0_sp=.001, b0_sp=.001,			# Hyperprior for tau.spatial
                       rd = .00001, # ratio of varSpike/varSlab
                       empiricalBayes = F,
                       dropZeros = F,
                       hCenter = F
){

  set.seed(seed) #random seed

  if(is.null(coord)){
    coord <- coord(msset)
  }



  if(dropZeros){
    res <- vector("list", length(feature))
    fInd <- 1
    seeds <- sample(size = length(feature), x=1:100000)

    for(f in feature){
      msset0 <- msset[f,]

      ######## drop zero pixels ##########
      for(s in levels(factor(techRep))){
        print(paste0("For feature", f, ", dropping ",
                     sum(spectra(msset0) == 0 & techRep == s),
                     " zero pixels from tissue ", s))
      }

      techRep0 <- techRep[c(spectra(msset0) != 0)]
      bioRep0 <- bioRep[c(spectra(msset0) != 0)]
      conditionOfInterest0 <- conditionOfInterest[c(spectra(msset0) != 0)]
      coord0 <- coord[c(spectra(msset0) != 0),]
      msset0 <- msset0[,c(spectra(msset0) != 0)]



      #### drop pixels with no neighbors ####
      W <- adj.grid(coord0, sample = factor(paste0(techRep0, conditionOfInterest0)),
                    type = type.neighbor,
                    radius = radius.neighbor,
                    max.dist = maxdist.neighbor)+0
      lonely <- which(rowSums(W) == 0)

      if(length(lonely) > 0){
        rm(W)
        print(paste0("dropping ", length(lonely), " pixel(s) with no neighbors: "))
        #print(coord0[lonely,])

        coord0 <- coord0[-lonely,]
        msset0 <- msset0[,-lonely]
        techRep0 <- techRep0[-lonely]
        bioRep0 <- bioRep0[-lonely]
        conditionOfInterest0 <- conditionOfInterest0[-lonely]

       if(any(table(techRep0) == 0)) warning("At least one tissue has no non-zero, non-island pixels")
        
        print("Remaining pixels per tissue:")
        print(table(techRep0))

        W <- adj.grid(coord0, sample = techRep0,
                      type = type.neighbor,
                      radius = radius.neighbor,
                      max.dist = maxdist.neighbor)+0
      }else{
        print("All pixels have at least one neighbor.")
      }

      ### Find any disconnected islands of pixels within samples (because spatial effects must be centered within islands) ####
      ### Code thanks to CARBayes ###
      W.list<- mat2listw(W)
      W.nb <- W.list$neighbours
      W.islands <- n.comp.nb(W.nb)
      islands <- W.islands$comp.id
      n.islands <- length(unique(islands))
      rm(W, W.list, W.nb, W.islands)

      if(length(levels(factor(techRep0))) > 1){
        if ( hCenter ) {
          print("Fitting model version: drop zeros, multiple samples, hiearchical centering")
          res[[fInd]] <- compareMSI_hc_zeros_multi(msset0,conditionOfInterest0,
                                          feature=1, nsim, burnin, trace,
                                          piPrior, seeds[fInd], logbase2, coord0,
                                          type.neighbor, radius.neighbor, maxdist.neighbor,
                                          spInit=NULL,
                                          bioRep0,
                                          techRep0,
                                          beta0, # Prior Mean for beta, only allow intercept
                                          prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                          precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                          a0_eps, b0_eps,     # Hyperprior for tau (1/eps.var)
                                          a0_bio, b0_bio,     # Hyperprior for taubio
                                          a0_tec, b0_tec,     # Hyperprior for tautec
                                          a0_sp, b0_sp,     # Hyperprior for tau.spatial
                                          rd, # ratio of varSpike/varSlab
                                          empiricalBayes,
                                          dropZeros = T
          )[[1]]
          res[[fInd]]$seed <- seeds[fInd]
        } else {
          print("Fitting model version: drop zeros, multiple samples, no hiearchical centering")
          res[[fInd]] <- compareMSI_multi(msset0,conditionOfInterest0,
                                          feature=1, nsim, burnin, trace,
                                          piPrior, seeds[fInd], logbase2, coord0,
                                          type.neighbor, radius.neighbor, maxdist.neighbor,
                                          spInit=NULL,
                                          bioRep0,
                                          techRep0,
                                          beta0, # Prior Mean for beta, only allow intercept
                                          prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                          precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                          a0_eps, b0_eps,     # Hyperprior for tau (1/eps.var)
                                          a0_bio, b0_bio,     # Hyperprior for taubio
                                          a0_tec, b0_tec,     # Hyperprior for tautec
                                          a0_sp, b0_sp,     # Hyperprior for tau.spatial
                                          rd, # ratio of varSpike/varSlab
                                          empiricalBayes,
                                          dropZeros = T
          )[[1]]
          res[[fInd]]$seed <- seeds[fInd]
        }
      }else{
        if ( hCenter ) warning("hierarchical centering not possible for n_tech = 1")
        print("Fitting model version: drop zeros, single sample, no hierarchical centering")
        res[[fInd]] <- compareMSI_zeros_single(msset0,conditionOfInterest0,
                                              feature=1, nsim, burnin, trace,
                                              piPrior, seeds[fInd], logbase2, coord0,
                                              type.neighbor, radius.neighbor, maxdist.neighbor,
                                              spInit=NULL,
                                              bioRep0,
                                              techRep0,
                                              beta0, # Prior Mean for beta, only allow intercept
                                              prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                              precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                              a0_eps, b0_eps,     # Hyperprior for tau (1/eps.var)
                                              a0_bio, b0_bio,     # Hyperprior for taubio
                                              a0_tec, b0_tec,     # Hyperprior for tautec
                                              a0_sp, b0_sp,     # Hyperprior for tau.spatial
                                              rd, # ratio of varSpike/varSlab
                                              empiricalBayes,
                                              dropZeros = T
        )[[1]]
        res[[fInd]]$seed <- seeds[fInd]
      }
      names(res)[fInd] <- paste0("Feature",f)
      fInd <- fInd + 1
    }
    return(res)

  }else{ #don't drop zeros

    techRep <- factor(techRep) #factor with different levels for each tissue (like "sample" before)
    n_tec <- length(levels(techRep)) #the number of distinct tissues
    nis_tec <- sapply(levels(techRep), function(x) sum(techRep == x)) #number of pixels in each tissue

    if(n_tec > 1){ 
      if ( hCenter ) { #do hiearchical centering for multi tissue experiments (only those without subsampling for now)
        print("Fitting model version: replace zeros with small value, multiple samples, hierarchical centering")
        return(compareMSI_hc_multi(msset,conditionOfInterest,
                                 feature, nsim, burnin, trace,
                                 piPrior, seed, logbase2, coord,
                                 type.neighbor, radius.neighbor, maxdist.neighbor,
                                 spInit,
                                 bioRep,
                                 techRep,
                                 beta0, # Prior Mean for beta, only allow intercept
                                 prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                 precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                 a0_eps, b0_eps,      # Hyperprior for tau (1/eps.var)
                                 a0_bio, b0_bio,      # Hyperprior for taubio
                                 a0_tec, b0_tec,      # Hyperprior for tautec
                                 a0_sp, b0_sp,      # Hyperprior for tau.spatial
                                 rd, # ratio of varSpike/varSlab
                                 empiricalBayes
        ))
      } else {
        print("Fitting model version: multiple samples, no hiearchical centering")
        return(compareMSI_multi(msset,conditionOfInterest,
                                 feature, nsim, burnin, trace,
                                 piPrior, seed, logbase2, coord,
                                 type.neighbor, radius.neighbor, maxdist.neighbor,
                                 spInit,
                                 bioRep,
                                 techRep,
                                 beta0, # Prior Mean for beta, only allow intercept
                                 prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                 precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                 a0_eps, b0_eps,      # Hyperprior for tau (1/eps.var)
                                 a0_bio, b0_bio,      # Hyperprior for taubio
                                 a0_tec, b0_tec,      # Hyperprior for tautec
                                 a0_sp, b0_sp,      # Hyperprior for tau.spatial
                                 rd, # ratio of varSpike/varSlab
                                 empiricalBayes
        ))
      }
    }else{ 
      if ( hCenter ) warning("hierarchical centering not possible for n_tech = 1")
       #don't do hiearchical centering if it's a single tissue experiment
        print("Fitting model version: single sample, no hierarchical centering")
        return(compareMSI_single(msset,conditionOfInterest,
                                 feature, nsim, burnin, trace,
                                 piPrior, seed, logbase2, coord,
                                 type.neighbor, radius.neighbor, maxdist.neighbor,
                                 spInit,
                                 bioRep,
                                 techRep,
                                 beta0, # Prior Mean for beta, only allow intercept
                                 prec0, # Prior Precision Matrix of beta (vague)  (only allow intercept)
                                 precAlpha0, #Prior Precision of slab (value of condition effect if it is not zero)
                                 a0_eps, b0_eps,      # Hyperprior for tau (1/eps.var)
                                 a0_bio, b0_bio,      # Hyperprior for taubio
                                 a0_tec, b0_tec,      # Hyperprior for tautec
                                 a0_sp, b0_sp,      # Hyperprior for tau.spatial
                                 rd, # ratio of varSpike/varSlab
                                 empiricalBayes
        ))
} #single or multi tissue experiment
  }#if we don't drop zeros and instead just add a small value before log transform

}#function



