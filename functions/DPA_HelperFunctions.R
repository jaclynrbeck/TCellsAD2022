# Functions used to run the DPA algorithm. "runDPA" will run all 3 steps in
# succession, or each of the 3 steps can be run separately by calling their
# individual functions. For testing use nIter = 100 or 1000, but for publication
# use at least nIter = 10,000. 

# Author: Jaclyn Beck
# Final script used for paper as of Sep 09, 2022

library(Seurat)
library(stringr)
library(scDC)
library(multcomp)
library(emmeans)

# Runs all 3 DPA steps in succession. Step 2 records the GLM summary (similar
# to an ANOVA) and Step 3 records post-hoc pairwise comparisons. 
# Required arguments:
# scRNA: Seurat object. Idents should be cluster identity. Cluster names should
#        be strings, NOT numbers. If cluster IDs are numbers, change them
#        to "Cluster <#>". 
# glmSummaryOutfile: full file path where the output of the pooled GLM 
#                    results will be saved. Should be a csv file. 
# pairwiseOutfile: full file path where the results for pairwise comparisons
#                  should go. Should be an xlsx file. 
# Optional arguments:
# step1intfile, step2intfile, step3intfile: full file paths where the 
#       intermediate results should be saved. Should be RDS files.
#       If NULL, intermediate result will not be saved. 
# nIter: number of bootstrap samples to generate. For the paper we used 10,000.
runDPA <- function( scRNA, glmSummaryOutfile, pairwiseOutfile, 
                    step1intfile = NULL, step2intfile = NULL, step3intfile = NULL, 
                    nIter = 10000 ) {
  res <- runDPAStep1(scRNA, step1intfile, nIter)
  fit <- runDPAStep2(res, glmSummaryOutfile, step2intfile)
  runDPAStep3(res, fit, pairwiseOutfile, step3intfile)
}


# Step 1: Generate bootstrapped samples.
# scRNA: Seurat object. Idents should be cluster identity. Cluster names should
#        be strings, NOT numbers. If cluster IDs are numbers, change them
#        to "Cluster <#>". 
# step1outfile: full file path where the intermediate result should be saved.
#               If NULL, intermediate result will not be saved. 
# nIter: number of bootstrap samples to generate.
# Returns "res", a data frame of bootstrapped samples
runDPAStep1 <- function( scRNA, step1outfile = NULL, nIter = 10000 ) {
  metadata <- data.frame(Sample = scRNA$orig.ident, 
                         Genotype = scRNA$genotype,
                         Cluster = Idents(scRNA))
  
  metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "5XFAD", "PS19", "PS-5X"))
  
  res = scDC_noClustering(cellTypes = metadata$Cluster,
                          subject = metadata$Sample,
                          calCI = TRUE,
                          calCI_method = c("multinom", "percentile"), 
                          nboot = nIter,
                          conf_level = 0.95,
                          ncores = 2,
                          verbose = TRUE)
  
  meta.info <- unique(metadata)
  res$results$genotype <- factor(str_replace(res$results$subject, "-[1|2]", ""),
                                 levels = levels(metadata$Genotype))
  res$info$genotype <- factor(str_replace(res$info$subject, "-[1|2]", ""),
                              levels = levels(metadata$Genotype))
  
  # Save intermediate result
  if (!is.null(step1outfile)) {
    saveRDS(res, step1outfile)
  }
  
  return(res)
}


# Step 2: Fit GLMs to each bootstrapped sample
# res: results data frame output from Step 1
# glmSummaryOutfile: full file path where the output of the pooled GLM 
#                    results will be saved. Should be a csv file. 
# step2outfile: full file path where the intermediate result should be saved.
#               If NULL, intermediate result will not be saved. 
# Returns "fit", a list containing individual GLM results and the pooled
#   GLM results.
runDPAStep2 <- function( res, glmSummaryOutfile, step2outfile = NULL) {
  minsample <- min(aggregate(value ~ subject, res$info, sum)$value)
  res$nstar <- round(res$thetastar * minsample) + 1
  
  fit <- fitGLM(res, res$info$genotype, subject_effect = FALSE, pairwise = FALSE, 
                fixed_only = TRUE, verbose = TRUE)
  #summary(fit$pool_res_fixed)
  subset(summary(fit$pool_res_fixed), p.value <= 0.05)
  
  # Save intermediate result
  if (!is.null(step2outfile)) {
    saveRDS(fit, step2outfile)
  }
  
  summ <- summary(fit$pool_res_fixed)
  summ$sig <- ' '
  summ$sig[summ$p.value <= 0.05] <- "*"
  write.csv(summ, glmSummaryOutfile, row.names = FALSE)
  
  return(fit)
}


# Step 3: Post-hoc individual comparisons
# res: results data frame output from Step 1
# fit: list containing "pool_res_fixed" (pooled GLM results) and "fit_fixed"
#      (all individual GLM results), output by Step 2.
# pairwiseOutfile: full file path where the results for pairwise comparisons
#                  should go. Should be an xlsx file. 
# step3outfile: full file path where the intermediate result should be saved.
#               If NULL, intermediate result will not be saved. 
runDPAStep3 <- function( res, fit, pairwiseOutfile, step3outfile = NULL ) {
  # Perform pairwise contrasts between genotypes in each cluster, for each 
  # individual GLM
  all_conts <- lapply(1:length(fit$fit_fixed), function(i) {
    glht(fit$fit_fixed[[i]], linfct=emm(pairwise~cond|cellTypes, adjust="tukey"))
  })
  
  # Get the coefficients from the above result
  all_ests <- lapply(1:length(all_conts), function(i) {
    lapply(1:length(all_conts[[i]]), function(x) {
      data <- summary(all_conts[[i]][[x]])$test$coefficients
    })
  })
  
  # Get the standard deviation from the above result
  all_ses <- lapply(1:length(all_conts), function(i) {
    lapply(1:length(all_conts[[i]]), function(x) {
      data <- summary(all_conts[[i]][[x]])$test$sigma
    })
  })
  
  # Save intermediate result
  if (!is.null(step3outfile)) {
    saveRDS(list(all_conts = all_conts, all_ests = all_ests, all_ses = all_ses),
            step3outfile)
  }
  
  # Turn coefficient and SD lists into one data frame
  clust_ests <- lapply(1:length(unique(res$info$cellTypes)), function(i) {
    data <- do.call(rbind, lapply(1:length(all_ests), function(a) {
      all_ests[[a]][[i]]
    }))
  })
  
  clust_ses <- lapply(1:length(unique(res$info$cellTypes)), function(i) {
    data <- do.call(rbind, lapply(1:length(all_ses), function(a) {
      all_ses[[a]][[i]]
    }))
  })
  
  # Pool results to get one significance value per pairwise comparison
  nmeans <- length(unique(res$info$genotype))
  nobs <- nrow(res$info)
  ncoefs <- nrow(summary(fit$pool_res_fixed))
  df <- nobs - ncoefs
  ncomp <- length(unique(res$info$cellTypes)) * choose(nmeans,2)
  
  sig <- lapply(1:length(clust_ests), function(C) {
    data <- do.call(rbind, lapply(1:ncol(clust_ests[[C]]), function(P) {
      pooled.s <- mice::pool.scalar(clust_ests[[C]][,P], clust_ses[[C]][,P]**2, n = nmeans)
      pooled.s <- t(unlist(pooled.s[4:length(pooled.s)])) 
      rownames(pooled.s) <- c(colnames(clust_ests[[C]])[P])
      pooled.s
    }))
    
    data <- as.data.frame(data)
    data$statistic <- data$qbar / sqrt(data$t)
    data$p.val <- 2*pt(-abs(data$statistic), df = df)
    data$p.val.adj <- p.adjust(data$p.val, method = "BH", n = ncomp)
    data$sig <- ' '
    data$sig[data$p.val.adj <= 0.05] <- "*"
    data <- cbind(rownames(data), data)
    colnames(data)[1] <- "comparison"
    data
  })
  
  names(sig) <- levels(res$info$cellTypes)
  sig
  
  writexl::write_xlsx(sig, path = pairwiseOutfile)
}

