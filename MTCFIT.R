
 MTCFIT <- function(mydata, id, trt, Y, typeY, vars_keep, myseed = 111,
                    ratio_k = c(1, 2, 3), caliper = c(0.1, 0.2, 0.3),
                    verbose = TRUE, alpha = 0.05, cut_length = 30) { 

  # ------------------ STOPPING CONDITIONS ------------------ #
  # stop if missing data
  if(all.equal(mydata[complete.cases(mydata),], mydata) == FALSE) stop("Cannot have missing data.") 
  
  # stop if user provides too many/not enough ratios
  if(length(ratio_k) < 3) {
    stop("Must have at least 3 choices for ratio_k.")
  }
  
  # stop if user provides too many/not enough calipers
  if(length(caliper) < 3) {
    stop("Must have at least 3 choices for caliper.")
  }
  
  # stop if wrong MatchIt version is being used
  suppressPackageStartupMessages(library(MatchIt, quietly = T))
    ActiveMatchItVers <- sessionInfo(package = "MatchIt")$otherPkgs[1]$MatchIt$Version
    if(ActiveMatchItVers != "3.0.2") stop(glue::glue("Must use MatchIt Version 3.0.2 and not version { ActiveMatchItVers }
                                                     The MatchIt archive can be accessed here: https://cran.r-project.org/web/packages/MatchIt/index.html."))
  # --------------------------------------------------------- #
  
  
  # -------------------- UPDATE OPTIONS --------------------- #
  old_options <- getOption("scipen")
  options(scipen = 9999)
  ## reset previous scipen setting on exit of MTCFIT
  on.exit(options(scipen = old_options))
  # --------------------------------------------------------- #
  
  
  # ---------------------------------------------------------#
  # Set the necessary libraries
  ############################################################
  
  suppressPackageStartupMessages(library(dplyr, quietly = T))
  suppressPackageStartupMessages(library(ggplot2, quietly = T))
  suppressPackageStartupMessages(library(ggpubr, quietly = T))
  suppressPackageStartupMessages(library(cowplot, quietly = T))
  suppressPackageStartupMessages(library(pROC, quietly = T))
  suppressPackageStartupMessages(library(MASS, quietly = T))
  suppressPackageStartupMessages(library(plyr, quietly = T)) 
  suppressPackageStartupMessages(library(doBy, quietly = T))
  suppressPackageStartupMessages(library(questionr, quietly = T))
  suppressPackageStartupMessages(library(buildmer, quietly = T))
  
  nvars <- length(vars_keep)
  
  # keep id, Y, and trt
  keep <- c(id, Y, trt)
  first3 <- mydata[keep]
  
  lastn <- mydata[vars_keep] 
  
  ## CREATE A LEGEND 
  regressors <- paste0("x", 1:nvars)
  lastn <- setNames(lastn, regressors)
  nomiss_data <- cbind(first3,lastn)

  set.seed(myseed)
  ##########################################################################
  # Create matrix of formulas for matching based on the number of variables
  ##########################################################################
  
  # construct permutation matrix of T/Fs of dimensions nvars^2 X nvars
  regMat <- expand.grid(rep(list(c(TRUE, FALSE)), nvars))
  # exclude intercept only model
  regMat <- regMat[apply(regMat, 1, sum) > 0, ]
  # order by number of regressors in the model
  regMat <- regMat[order(apply(regMat, 1, sum)), ]
  # construct list containing all possible formulas (save for intercept only, "trt ~ 1")
  allMatchingList <- apply(regMat, 1, function(x) arsenal::formulize(trt, x = regressors[x]))

  crit_value <- qnorm(1 - alpha/2)
  conf_level <- 1-alpha
  conf_level_label <- paste0(round( 100 * conf_level, 2), "%")
  
  #################################################################################################################
  #prep for computational model - match on any significantly associated variable with the treatment in the dataset
  #################################################################################################################
  
  # drop id and y from initial dataset and run model to see what is associated with trt
  dropthese <- names(mydata) %in% c(id, Y) # converted to variable name
  newmydt <- mydata[!dropthese]
  
  # Fit the model
  allpm <- glm(arsenal::formulize(y = trt, x = "."), data = newmydt, family="binomial") %>%
    stepAIC(trace = FALSE)
  
  # Summarize the final selected model
  if(verbose) summary(allpm) 
 
  
  #################################################################################################################
  #Computational model - evaluate a match using all variables in dataset
  #################################################################################################################
  
  # prep for computational model - match on any significantly associated variable with the treatment in the dataset
  # drop id and y from initial dataset and run model to see what's associated with trt
  
  # Fit the model with stepwise
  step_trace <- 1
  if(verbose == FALSE) step_trace <- FALSE
  step1 <- step(glm(arsenal::formulize(trt, "."), data = newmydt, family = "binomial"), 
                direction = 'both', trace = step_trace)
  
  match_formula <- formula(step1)
  match_formula_rhs <- match_formula[[3]]
 
  results_comp <- matrix(nrow = 0, ncol = 11)
  
  ########################################################
  est_effect <- function(fn, method, .alpha = alpha, 
                         .match_formula = match_formula, 
                         .conf_level = conf_level) {
    
    crit_val <- qnorm(1 - .alpha/2) 
    
    mtc <- match.data(fn)
    
    # Estimating the treatment effect + xx% CI
    e <- glm(arsenal::formulize(Y, trt), data = mtc, family = typeY)
    AIC <- AIC(e)
    
    eff <- summary(e)$coefficients[2,1]
    sd <- summary(e)$coefficients[2,2]

    if (typeY == "binomial") {
      eftr = exp(eff)
      ci_low <- exp(eff - crit_val * sd) 
      ci_high <- exp(eff + crit_val * sd)
    } else {
      eftr <- eff * 1
      ci_low <- eff - crit_val * sd
      ci_high <- eff + crit_val * sd
    }
  
    # Estimating the AUC for the matching model
    tde <- glm(.match_formula, data = mtc, family=typeY) 
    pr <- data.frame(predict(tde, type = 'response'))
    a_e <-cbind(mtc, pr)
    roc_obj <- suppressMessages(roc(a_e[[trt]], a_e$pr))
    area <- auc(roc_obj)
    ci_auc_lower <- as.numeric(ci.auc(area, conf.level = conf_level)[1])
    ci_auc_upper <- as.numeric(ci.auc(area, conf.level = conf_level)[3])
    
    method <- c(method)
    combine <- data.frame(method, eftr, sd, ci_low, ci_high, k, cal, AIC, 
                          area, ci_auc_lower, ci_auc_upper)
    
    combine
  }
  ########################################################
  
  for (h in seq_along(ratio_k))      {
    
    k <- ratio_k[h]
    
    for (l in seq_along(caliper))   {
      
      cal <- caliper[l]
    
      # for computational;
      # Nearest Neighbor with replacement, and PS distance, target estimand ATT
      set.seed(myseed)
      fn1_ <- matchit(match_formula, data = mydata, replace = T,
                    method = "nearest", caliper = cal, ratio = k)
      
      # Glm  with replacement and Mahalanobis distance matching,target estimand ATT
      set.seed(myseed)
      fn2_ <- matchit(match_formula, data = mydata, replace = T, mahvars = ~ match_formula_rhs,
                      distance = "mahalanobis", caliper = cal, ratio = k)
  
      combine3 <- est_effect(fn1_, "Nearest Neighbor")
      combine4 <- est_effect(fn2_, "Mahalanobis")
      
      allcomb_ <- rbind (combine3, combine4)
      results_comp <- rbind(results_comp, allcomb_) 
    
    }
  }
  
  results_comp$PredsIn <- c(match_formula_rhs)
  
  ##########################################################################
  # Matching pn two methods with all combinations of method, caliper, ratio - PI selected X's
  ##########################################################################
  #s et caliper and ratio for matching
  
  results <- matrix(nrow = 0, ncol = 17) ## ncol = 16
 
  # unadjusted effect + xx% CI
  ujs <- glm(arsenal::formulize(Y, trt), data = nomiss_data, family = typeY)
  
  if(verbose) summary(ujs)
  
  unadj_eff <- summary(ujs)$coefficients[2, 1]
  unadj_sd <- summary(ujs)$coefficients[2, 2]
  
  if (typeY == "binomial") {
    unadj_efftr = exp(unadj_eff)
    unadj_low <- exp(unadj_eff - crit_value * unadj_sd)
    unadj_high <- exp(unadj_eff + crit_value * unadj_sd)
    } else {
    unadj_efftr <- unadj_eff * 1
    unadj_low <- unadj_eff - crit_value * unadj_sd
    unadj_high <- unadj_eff + crit_value * unadj_sd
    }
  
  ########################################################
  est_effect2 <- function(fn, method, mtv, 
                          .alpha = alpha, .conf_level = conf_level) { 
    
    mtc <- match.data(fn)
    
    crit_value <- qnorm(1-.alpha/2)
   
    # Estimating the treatment effect + 95% CI
    e <- glm(arsenal::formulize(Y, trt), data = mtc, family = typeY)
    AIC <- AIC(e)
    
    eff <- summary(e)$coefficients[2, 1]
    sd <- summary(e)$coefficients[2, 2]
    p <- summary(e)$coefficients[2, 4]
    
    if (typeY == "binomial") {
      eftr = exp(eff)
      ci_low <- exp(eff - crit_value * sd)
      ci_high <- exp(eff + crit_value * sd)
    } else {
      eftr <- eff * 1
      ci_low <- eff - crit_value * sd
      ci_high <- eff + crit_value * sd
    }
    
    # Estimating the AUC for the matching model
    tde <- glm(mtv, data = mtc, family = "binomial")
    pr <- data.frame(predict(tde, type = 'response'))
    a_e <- cbind(mtc, pr)
    roc_obj <- suppressMessages(roc(a_e[[trt]], a_e$pr))
    area <- auc(roc_obj)
    ci_auc_lower <- as.numeric(ci.auc(area, conf.level = conf_level)[1])
    ci_auc_upper <- as.numeric(ci.auc(area, conf.level = conf_level)[3])
    
    method <- c(method)
    PredsIn <- deparse(mtv)
    combine <- data.frame(method, eftr, sd, p, ci_low, ci_high, k, cal, 
                          unadj_efftr, unadj_sd, unadj_low, unadj_high, 
                          PredsIn, AIC, area, ci_auc_lower, ci_auc_upper)
    combine
  }
  ########################################################
  
  for(h in seq_along(ratio_k)) {
    
    k <- ratio_k[h]
    
    for(l in seq_along(caliper)) {
      
      cal <- caliper[l]
      
      for(v in seq_along(allMatchingList)) {
        
        mdl <- allMatchingList[[v]]
        
        mah <- mdl[[3]]
        
       
        # Nearest Neighbor with replacement, and PS distance, target estimand ATT
        set.seed(myseed)
        fn1 <- matchit(mdl, data = nomiss_data, replace = T,
                     method = "nearest", caliper = cal, ratio = k)
        
        # GLM  with replacement and Mahalanobis distance matching,target estimand ATT
        set.seed(myseed)
        fn2 <- matchit(mdl, data = nomiss_data, replace = T, mahvars = ~ mah,
                      distance = "mahalanobis", caliper = cal, ratio = k) # glm
    
        combine1 <- est_effect2(fn1, "Nearest Neighbor", mdl)
        combine2 <- est_effect2(fn2, "Mahalanobis", mdl)
        
        allcomb <- rbind(combine1, combine2)
        
        results <- rbind(results, allcomb)  
        
      }
    }
  }
  
  # custom theme for plots
  theme_set(theme_cowplot())
  
  mytheme = list(
    theme_classic() +
      theme(panel.background = element_blank(), strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(), strip.text = element_text(face="bold", size=7),
            axis.text=element_text(face="bold"), axis.title = element_text(face="bold"),
            plot.title = element_text(face = "bold", hjust = 0.5,size=9)))
  
  # colors
  best_rand_lty <- "longdash" # or lty = 5
  vert_bar_color <- "#e72b56" # rm "...FF"
  ref_line_color <- "black"
  NN_color <- "#a21fa8"
  MAH_color <- "#d4840d"
  methods_colors <- c(NN_color, MAH_color)
  

  ##############################################################################
  #GRAPH 1 - estimated effect + xx%CI for match on all variables by method+K+cal 
  #with unadjusted effect +95%error band
  ##############################################################################
  
  ttt <- deparse(allMatchingList[[which(apply(regMat, 1, sum) == nvars)]])
  results_gr1 <- results %>% filter(PredsIn == ttt)
  
  diffl <- unadj_high - unadj_low
  
  # label for plots
  ttt_x_label <- sub(x = ttt, pattern = trt, replacement = "trt")
  
  uplab <- unadj_efftr + .02
 
  # position for label in graph
  if (typeY == "binomial") { 
    uplab_up <- 1.02
  } else {
    uplab_up <- .02
  }
  
  
  
  # combine k and cal for graph
  results_gr1$kc <- as.factor(paste(results_gr1$k, results_gr1$cal, sep = " & "))
  
  labeleff <- ifelse(typeY == "binomial", 
                     paste0("Odds Ratio + ", conf_level_label, " CI"), 
                     "Coefficient")
  
  label_eff_line <- ifelse(typeY == "binomial", "Effect = 1", "Effect = 0")
  
  eff_line_val <- ifelse(typeY == "binomial", 1, 0)
  
  combon <- nvars * (nvars - 1)/4
  
  
  myplot1 <- ggplot(data = results_gr1, aes(x = kc, y = eftr, color = method)) +
    geom_point(aes(kc, eftr), size = 1, position = position_dodge(.2)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high, linetype = method), 
                  width = .5, position = position_dodge(.2), size = .5) +
    scale_x_discrete(name = "ratio & caliper", guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = labeleff) + 
    geom_hline(aes(yintercept = unadj_efftr), results_gr1) + 
    annotate(geom = "text", x = combon, y = uplab, label = "Unmatched effect", size = 2.3) + 
    geom_hline(aes(yintercept = unadj_low), results_gr1, linetype = "dashed") +
    geom_hline(aes(yintercept = unadj_high), results_gr1, linetype = "dashed") + 
    geom_hline(aes(yintercept = eff_line_val), linetype = "dotdash") + # 1
    annotate(geom = "text", x = 2.5, y = uplab_up, label = label_eff_line, size = 2.3) + # Effect=1
    mytheme  +
    scale_color_manual(values = c(NN_color, MAH_color))
  
  ##################################################################
  #GRAPH 2 - Median effect + xx%CI from all methods+cal+k 
  #with effect from match on all variables + xx%error band
  #in this step the best matching formula is selected
  ##################################################################
  #pull median OR and CI for each model;
  S_med_OR <- data.frame(summaryBy(eftr ~ PredsIn + method, data = results, 
                                 FUN = list(median))) # median of n(ratio) X n()
  S_med_LL <- data.frame(summaryBy(ci_low ~ PredsIn + method, data = results, 
                                 FUN = list(min)))
  S_med_LL <- subset(S_med_LL, select = -c(PredsIn, method))
  S_med_UL <- data.frame(summaryBy(ci_high ~ PredsIn + method, data = results, 
                                 FUN = list(max)))
  S_med_UL <- subset(S_med_UL, select = -c(PredsIn, method))
  
  S_med_all_ <- cbind(S_med_OR, S_med_LL, S_med_UL)
  S_med_all<- S_med_all_ %>% arrange(desc(PredsIn))
  
  S_med_all$ID <- seq.int(nrow(S_med_all))
  
  
  #prep for graph using match on all predictor as referent
  S_med_red <- subset(S_med_all, S_med_all$ID != c(1,2))
  S_med_lab1 <- subset(S_med_all, S_med_all$ID == 1)
  S_med_lab2 <- subset(S_med_all, S_med_all$ID == 2)
  
  
  #3 3=MEDIAN, 4=LL, 5=UL
  colnames(S_med_lab1) <- paste0("NN", 1:ncol(S_med_lab1))

  colnames(S_med_lab2) <- paste0("MA", 1:ncol(S_med_lab2))

  # position for label in graph
  if (S_med_lab1$NN3 > S_med_lab2$MA3) {
    uplab1 <- S_med_lab1$NN3 + .02
    uplab2 <- S_med_lab2$MA3 - .02
  } else {
    uplab1 <- S_med_lab1$NN3 - .02
    uplab2 <- S_med_lab2$MA3 + .02
  }
  
  
  # label PredsIn
  S_med_red$PredsIn <- sub(x = S_med_red$PredsIn, pattern = trt, replacement = "trt")
  
  myplot2 <- ggplot(data = S_med_red, aes(x = PredsIn, y = eftr.median, color = method)) +
    geom_point(aes(PredsIn, eftr.median), size = 1, position = position_dodge(.2)) +
    geom_errorbar(aes(ymin = ci_low.min, ymax = ci_high.max, linetype = method), width = .5, 
                  position = position_dodge(.2), size = .5) +
    scale_x_discrete(name = "Matching Method and Equation", 
                     guide = guide_axis(angle = 45)) + 
    scale_y_continuous(name = labeleff) + 
    geom_hline(aes(yintercept = NN3), S_med_lab1, colour = NN_color) +
    geom_hline(aes(yintercept = NN4), S_med_lab1, linetype = "dashed" , 
               colour = NN_color) +
    geom_hline(aes(yintercept = NN5), S_med_lab1, linetype = "dashed", 
               colour = NN_color) +
    annotate(geom = "text", x = combon, y = uplab1, label = ttt_x_label, size = 2,  
             colour = NN_color) +
    geom_hline(aes(yintercept = MA3), S_med_lab2, colour = MAH_color) +
    geom_hline(aes(yintercept = MA4), S_med_lab2,  linetype = "dashed", 
               colour = MAH_color) +
    geom_hline(aes(yintercept = MA5), S_med_lab2,  linetype = "dashed", 
               colour = MAH_color) +
    annotate(geom = "text", x = combon, y = uplab2, label = ttt_x_label, size = 2.3, 
             colour = MAH_color) +
    geom_hline(aes(yintercept = eff_line_val), linetype = "dotdash") + # 1
    annotate(geom = "text", x = 2.5, y = uplab_up, label = label_eff_line, size = 2.3) + # "Effect=1"
    mytheme +
    scale_color_manual(values = c(NN_color, MAH_color))

  # calculate number of ORs >1, =1, <1;
  if(typeY == "binomial") {
    results_forp <- results %>% mutate(
      group_var = case_when(
        p < alpha & eftr < 1 ~ paste0("OR<1 (p<", alpha, ")"),
        p < alpha & eftr > 1 ~ paste0("OR>1 (p<", alpha, ")"),
        TRUE ~ paste0("OR=1 (p>", alpha, ")") 
        )
    )
  } else {
    results_forp <- results %>% mutate(
      group_var = case_when(
        p < alpha ~ paste0("Coeff not = 0 (p<", alpha, ")"),
        TRUE ~ paste0("Coeff = 0 (p>", alpha, ")") 
      )
    )
  }
    
  
  ########################################################
  #######################################################
  
  tab1NAout <- subset(results_forp, method == "Nearest Neighbor")
  tab1MAout <- subset(results_forp, method == "Mahalanobis")
  tab1NA <- as.data.frame(summary(arsenal::tableby(group_var ~ PredsIn, 
                                                   cat.stats = "countrowpct", 
                                                   total = FALSE,
                                                   test = FALSE,
                                                   data = tab1NAout),
                                  text = "html"))
  
  names(tab1NA)[2] <- tab1MAout$group_var[1] 
  
  # keep everything in the var names up to and including the first ")"
  old_tab1NA_names <- names(tab1NA)
  new_tab1NA_names <- old_tab1NA_names %>% stringr::str_extract(., "^.*?\\)")
  names(tab1NA) <- new_tab1NA_names
    
  tab1NA <- tab1NA[-c(1), ]
  type <- stringr::str_remove_all(tab1NA[, 1], "[&nbsp;]")
  method <- "Nearest Neighbor"
  
  tab1NA <- tab1NA[-c(1)]
  
  if(ncol(tab1NA) == 1) {
    tab1NA[, 1] <- tab1NA[, 1] %>% stringr::str_extract("(?<=\\().+?(?=\\))")
  } else {
    tab1NA <- tab1NA %>% 
      apply(., 2, stringr::str_extract, pattern = "(?<=\\().+?(?=\\))")
  }
  
  tab1NA <- cbind(method, type, tab1NA)

  
  tab1MA <- as.data.frame(summary(arsenal::tableby(group_var ~ PredsIn, 
                                                   cat.stats = "countrowpct", 
                                                   total = FALSE,
                                                   test = FALSE,
                                                   data = tab1MAout),
                                  text = "html"))
  
  # keep everything in var names up until and including the first  ")"
  old_tab1MA_names <- names(tab1MA)
  new_tab1MA_names <- old_tab1MA_names %>% stringr::str_extract(., "^.*?\\)")
  names(tab1MA) <- new_tab1MA_names
  
  tab1MA <- tab1MA[-c(1), ]
  type <- stringr::str_remove_all(tab1MA[, 1], "[&nbsp;]")
  method <- "Mahalanobis"

  tab1MA <- tab1MA[-c(1)]
  
  if(ncol(tab1MA) == 1) {
    tab1MA[, 1] <- tab1MA[, 1] %>% stringr::str_extract("(?<=\\().+?(?=\\))")
  } else {
    tab1MA <- tab1MA %>% 
      apply(., 2, stringr::str_extract, pattern = "(?<=\\().+?(?=\\))")
  }
  
  tab1MA <- cbind(method, type, tab1MA)
  
  ######################################################################################
  #GRAPH 3 - Best matching formula effect + each of 5 random variables from dataset
  #with best matching formula effect + xx%error band
  #in this step compare "best matching formula" to best of "best matching formula+random"
  #######################################################################################
  #select best method based on:
  #1 - if OR+CI are defined
  #2 - largest AUC
  #3 - smallest SD 
  #4 - smallest num of preds
  #5 - most important predictor
  
  results$eftr <- round(results$eftr, 2)
  results$ci_low <- round(results$ci_low, 2)
  results$ci_high <- round(results$ci_high, 2)
  results$AIC <- round(results$AIC, 2)
  results$area <- round(results$area, 2)
  results$areall <- round(results$ci_auc_lower, 2)
  results$areaul <- round(results$ci_auc_upper, 2)
  
  #verify and clean from indetermined ORs or CI
  # 1
  results_2 <- results %>%  
    filter_all(all_vars(is.finite(.))) 
  
  # 3-5
  results_sort <- results_2 %>% arrange(desc(area), sd, PredsIn)
  thebest <- results_sort[1, ]
  
  # subset original dataset for random selection of 5 addditional matching variables 
  # exlcude id,y,trt and original matching variables
  
  df <- mydata[, !(names(mydata) %in% vars_keep)]
  dataforZ <- df[, !(names(df) %in% keep)]
  
  # randomly select 5 additional matching variables
  random_sel <- ifelse(length(dataforZ) < 5, length(dataforZ), 5)
  
  # stop if additional variable criteria is not met.
  if(random_sel < 1) stop("Must have at least one additional variable for random selection.")
  
  Z_kp <- sample(colnames(dataforZ), random_sel)
  newpool <- mydata[Z_kp]
  
  # rename predictors
  names(newpool) <- paste0("z", 1:random_sel)
  
  # attached to set with X's and keep the clean original version (no missing)
  mydfz <- cbind(nomiss_data, newpool)
  
  #keep parameters for last function 
  tbfml <- as.character(thebest$PredsIn)
  tbmtd <- as.character(thebest$method)
  trasnf <- as.formula(tbfml)
  tbcal <- thebest$cal
  tbk <- thebest$k
  tbef <- thebest$eftr
  tbll <- thebest$ci_low
  tbul <- thebest$ci_high
  tbarea <- thebest$area
  tbareall <- thebest$areall
  tbareaul <- thebest$areaul
  tbaic <- thebest$AIC
  tbsd <- thebest$sd
  
  #implement formula with Z's and re-run match
  lastfun <- function (Z, tbcal, tbk, .alpha = alpha, 
                       .conf_level = conf_level) {
    
    frmlnew <- add.terms(trasnf, Z)
    
    if (tbmtd == "Nearest Neighbor") {
      rptfn <- matchit(frmlnew, data = mydfz, replace = TRUE, 
                       method = "nearest", caliper = tbcal, ratio = tbk)
      } else {
        rptfn <- matchit(frmlnew, data = mydfz, replace = TRUE, 
                          distance = "mahalanobis", caliper = tbcal, ratio = tbk)
      }
    
    mtc_dt <- match.data(rptfn)
    crit_value <- qnorm(1 - .alpha/2)
    
    # Estimating the treatment effect + xx% CI
    e <- glm(arsenal::formulize(Y, trt), data = mtc_dt, family = typeY) # does all of this need to be piped in as well?
    AIC <- AIC(e)
    
    eff <- summary(e)$coefficients[2, 1]
    sd <- summary(e)$coefficients[2, 2]
    
    if (typeY == "binomial") {
      eftr = round(exp(eff), 2)
      ci_low <- round(exp(eff - crit_value * sd), 2)
      ci_high <- round(exp(eff + crit_value * sd), 2)
    } else {
      eftr = round(eff * 1, 2) 
      ci_low <- round(eff - crit_value * sd, 2)
      ci_high <- round(eff + crit_value * sd, 2)
    }
  
    #Estimaing the AUC for the matching model
    tde <- glm(frmlnew, data = mydfz, family = "binomial")
    pr <- data.frame(predict(tde, type = 'response'))
    a_e <- cbind(mydfz, pr)
    roc_obj <- suppressMessages(roc(a_e[[trt]], a_e$pr))
    area <- auc(roc_obj)
    ci_auc_lower <- as.numeric(ci.auc(area, conf.level = .conf_level)[1])
    ci_auc_upper <- as.numeric(ci.auc(area, conf.level = .conf_level)[3])
    
    method <- c(tbmtd)
    prds <- deparse(frmlnew)
    totvr <- data.frame(method, eftr, sd, ci_low, ci_high, tbk, tbcal, prds, AIC, 
                        Z, tbef, tbsd, tbll, tbul, tbaic, area, ci_auc_lower, 
                        ci_auc_upper)
    totvr
  }

 
  with_Z1 <- lastfun("z1", tbcal, tbk)
  with_Z2 <- lastfun("z2", tbcal, tbk)
  with_Z3 <- lastfun("z3", tbcal, tbk)
  with_Z4 <- lastfun("z4", tbcal, tbk)
  with_Z5 <- lastfun("z5", tbcal, tbk)
  
  finalZ <- rbind(with_Z1, with_Z2, with_Z3, with_Z4, with_Z5)
  
  #position for label in graph
  upzip1 <- finalZ$tbef + .03
  
  diffz2 <- tbul - tbll
  
  clr <- ifelse(tbmtd == "Mahalanobis", MAH_color, NN_color)
  
 
  # label for plots
  tbfml_x_label <- sub(x = tbfml, pattern = trt, replacement = "trt")
  
  myplot3 <- ggplot(data = finalZ, aes(x = Z, y = eftr)) +
    geom_point(aes(Z, eftr),  size = 1, position = position_dodge(.2)) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high, linetype = method), width = .5, position = position_dodge(.2), size = .5) +
    scale_x_discrete(name = "Best Match with Randomly Selected Covariate", guide = guide_axis(angle = 45)) + 
    scale_y_continuous(name = labeleff) + 
    geom_hline(aes(yintercept = tbef), finalZ, color = clr) + 
    annotate(geom = "text", x = combon, y = upzip1, label = tbfml_x_label, size = 2.5, color = clr) + 
    geom_hline(aes(yintercept = tbll), finalZ,  linetype = "dashed", colour= clr) +
    geom_hline(aes(yintercept = tbul), finalZ,  linetype = "dashed", colour= clr) +
    geom_hline(aes(yintercept = eff_line_val), linetype = "dotdash") + # 1
    annotate(geom = "text", x = 2.5, y = uplab_up, label = label_eff_line, size = 2.3) + 
    mytheme +
    theme(legend.title = element_blank()) + theme(legend.position = "none")
  
  ###########################################
  # identify best final model with random var
  ###########################################
  
  #verify and clean from indetermined ORs or CI
  #1
  finalZ_2 <- finalZ %>% 
    filter_all(all_vars(is.finite(.))) %>%
    mutate(prds = sub(x = prds, pattern = trt, replacement = "trt"))
  
  #3-5
  finalZ_3 <- finalZ_2 %>% arrange(desc(area), sd, prds)
  finalZfZ <- finalZ_3[1, ]
  
  Zsortf <- finalZfZ[order(finalZfZ$AIC, finalZfZ$sd), ]
  thebestwithz <- Zsortf[1, ]
  
  prdsfnz <- as.character(thebestwithz$prds)
  tbeffnz <- round(thebestwithz$eftr, 2)
  tbllfnz <- round(thebestwithz$ci_low, 2)
  tbulfnz <- round(thebestwithz$ci_high, 2)
  tbaicfnz <- round(thebestwithz$AIC, 2)
  tbareafnz <- round(thebestwithz$area, 2)
  tbareallfnz <- round(thebestwithz$ci_auc_lower, 2)
  tbareaulfnz <- round(thebestwithz$ci_auc_upper, 2)
  
  ###########################################
  # dentify the best computational
  ###########################################
  # verify and clean from undetermined ORs or CI
  results_comp$eftr <- round(results_comp$eftr, 2)
  results_comp$ci_low <- round(results_comp$ci_low, 2)
  results_comp$ci_high <- round(results_comp$ci_high, 2)
  results_comp$AIC <- round(results_comp$AIC, 2)
  results_comp$area <- round(results_comp$area, 2)
  results_comp$areall <- round(results_comp$ci_auc_lower, 2)
  results_comp$areaul <- round(results_comp$ci_auc_upper, 2)
  
  # 3-5
  results_comp2so <- results_comp %>% arrange(desc(area), sd)
  comp <- results_comp2so[1, ]
  
  compfml <- as.character(comp$PredsIn)
  compmtd <- as.character(comp$method)
  compcal <- comp$cal
  compk <- comp$k
  compef <- comp$eftr
  compll <- comp$ci_low
  compul <- comp$ci_high
  comparea <- comp$area
  compareall <- comp$areall
  compareaul <- comp$areaul
  compaic <- comp$AIC
  compsd <- comp$sd
  
  
  ###########################################
  # OUTPUT
  ##########################################
  tlabx1 <- paste(vars_keep, collapse = ",  ")
  tlabx2 <- paste(Z_kp, collapse = ",  ")
  Zs <- paste0("z", 1:5)
  
  new_table_body_df <- data.frame(`Match Type` = c("Best", "Random", "Computational"),
                                  `Method (cal, ratio)` = c(
                                    paste0(tbmtd," (", tbcal, ", ", tbk, ")"), 
                                    paste0(as.character(thebestwithz$method), " (",  # this best method, ratio, cal same as "best"
                                           thebestwithz$tbcal, ", ", thebestwithz$tbk, ")"),
                                    paste0(compmtd, " (", compcal, ", ", compk, ")")
                                    ),
                                  Model = c(stringr::str_wrap(tbfml, 20),
                                            stringr::str_wrap(prdsfnz, 20),
                                            stringr::str_wrap(paste0(trt, " ~ ", paste0(compfml, collapse = " + ")), 20)
                                            ),
                                  AUC = c(paste0(tbarea, " (", tbareall, ", ", tbareaul, ")"),
                                          paste0(tbareafnz, " (", tbareallfnz, ", ", tbareaulfnz, ")"),
                                          paste0(comparea, " (", compareall, ", ", compareaul, ")")
                                          ),
                                  Estimate = c(paste0(tbef," (", tbll, ", ", tbul, ")"),
                                               paste0(tbeffnz," (", tbllfnz, ", ", tbulfnz, ")"),
                                               paste0(compef, " (", compll, ", ", compul, ")")
                                               ),
                                  AIC = c(tbaic, tbaicfnz, compaic),
                                  check.names = FALSE
                                  ) 
  
  names(new_table_body_df)[c(4, 5)] <- c(paste0("AUC (", conf_level_label, " CI)"),
                                         paste0("Estimate (", conf_level_label, " CI)"))
                                         
  
  new_table_body <- new_table_body_df %>%
    ggtexttable(theme = ttheme(tbody.style = tbody_style(size = 10, 
                                                         hjust = 0, 
                                                         x = 0.1)),
                rows = NULL) 

  length_of_formula <- function(x) {
    nterms <- length(x)
    sumchar <- 0
    cumchar <- vector(mode = "numeric")
    for(fterm in seq_along(x)) {
      sumchar <- nchar(x[fterm]) + sumchar
      cumchar[fterm] <- sumchar
    }
    cumchar
  }
  
  cut_string_here <- function(x, .cut_length = cut_length) {
    which.min(abs(x - .cut_length))
  }
  
  paste_split_formula <- function(x, cut_location) {
    if(cut_location != length(x)) {
      paste0(
        paste0(x[1:cut_location], collapse = ", "), 
        "\n\t\t", 
        paste0(x[(cut_location + 1):length(x)], collapse = ", ")
      )
    } else {
      paste(x, collapse = ", ")
      }
    }
  
  # split initial matching
  split_vars_keep <- cut_string_here(length_of_formula(vars_keep))
  initial_matching_vec <- stringr::str_c(regressors, rep(" = ", length(regressors)), vars_keep)
  initial_matching_text <- paste_split_formula(initial_matching_vec, split_vars_keep)
  
  # split randomly selected
  split_Z_kp <- cut_string_here(length_of_formula(Z_kp))
  random_matching_vec <- stringr::str_c(Zs, rep(" = ", length(Zs)), Z_kp)
  random_matching_text <- paste_split_formula(random_matching_vec, split_Z_kp)

  new_table_subtitle_text <- glue::glue(
  "
  Outcome:                                    { Y }
  Treatment:                                   { trt }
  Initial Matching Covariates:           { initial_matching_text }
  Randomly Selected Covariates:    { random_matching_text }
  Computational Model Covariates:  All covairates in { substitute(mydata) }.
  "
  )

  bottom_right_panel <- new_table_body %>%
    tab_add_title(text = new_table_subtitle_text, face = "plain", size = 10) # make smaller size to account for additoinal width
  
  left_panel <- ggarrange(myplot1, myplot2 + font("x.text", size = 9),
                          common.legend = TRUE, legend = "bottom",
                          ncol = 1, nrow = 2,
                          heights = c(1, 1))

  right_panel <- ggarrange(myplot3, bottom_right_panel,
                           ncol = 1, nrow = 2,
                           heights = c(1, 1),
                           widths = c(1, 0.75))
  
  final_plot <- ggarrange(left_panel, right_panel, ncol = 2, nrow = 1)
  supp_table <- rbind(tab1NA, tab1MA)
  
  out <- list(plot_table = new_table_body_df, supplemental_table = supp_table, 
              plot = final_plot)
  
  print(out)
}
