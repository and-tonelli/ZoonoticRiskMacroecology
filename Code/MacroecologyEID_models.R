# ------------------------------------------------------------------------------
# 'Advancing the use of macroecological approaches for the prediction of zoonotic disease risk'
# Corresponding authors: moreno.dimarco@uniroma1.it; andrea.tonelli@uniroma1.it

###load packages
library(randomForestSRC); library(caTools); library(tidyr); library(dplyr); library(ggplot2); 
library(gridExtra); library(forcats); library(corrplot); library(grid)

# Load main dataframe with environmental and macroecological variables
load("Data/MacroecologyEID_df.RData")
df <- na.omit(df) # omit municipalities for which we have no data

#reordering dataframe 
d_cutleish <- df[,c("Fold", "Human_footprint","Population","Cropland", "Forest","Pasture","Annual_tmp","Annual_pre",
                    "Annual_wet","BII","Host_rich", "EVI_varian","Forest_los", "Fun_rich", "Output"
)]

d_cutleish$Output <- as.factor(as.character(d_cutleish$Output))
d_cutleish$Output <- fct_relevel(d_cutleish$Output, "1")


# Correlation ####

d_corr <- d_cutleish[,c("Human_footprint","Population","Cropland","Forest",
                        "Pasture","Annual_tmp","Annual_pre","Annual_wet",
                        "BII","Host_rich", "EVI_varian", "Forest_los", 
                        "Fun_rich"
)]

colnames(d_corr) <- c("Human footprint", "log(Population)", "Cropland", "Forest (cover)", 
                      "Pasture", "Temperature","Precipitation", "Wet days", 
                      "BII", "Host richness", "EVI", "Forest loss", 
                      "Functional richness"
)


Pearson_cor_matrix <- cor(d_corr, method = "pearson")
Spearman_cor_matrix <- cor(d_corr, method = "spearman")

# plot correlation
corrplot(
  Pearson_cor_matrix, method = 'circle', order = 'original',  type = 'lower', diag = FALSE,
  tl.cex = 1.2, tl.col = "black", cl.cex = 0.9, addCoef.col = "black", number.cex = 1,
  col = colorRampPalette(c("#f3722c", "#f9844a", "#edf2f4", "#81c3d7", "#3a7ca5"))(200)
)

corrplot(
  Spearman_cor_matrix, method = 'circle', order = 'original',  type = 'lower', diag = FALSE,
  tl.cex = 1.2, tl.col = "black", cl.cex = 0.9, addCoef.col = "black", number.cex = 1,
  col = colorRampPalette(c("#f3722c", "#f9844a", "#edf2f4", "#81c3d7", "#3a7ca5"))(200)
)

t1 <- Sys.time()

#### BASE MODEL ####

#creating formula
n <- setdiff(names(d_cutleish), c("Fold", "BII","Host_rich", "EVI_varian", "Forest_los", "Fun_rich"))

cl_formula <- as.formula(paste("Output~", paste(n[-length(n)], collapse = " + ")))


##------------------------------------------------##
##   for loop for cutaneous leishmaniasis model   ##
##------------------------------------------------##
oob_out_cv <- data.frame(matrix(vector(),0, 6, dimnames=list(c(), c("auc", "sens", "spec", "TSS", "Type",  "spatial_fold"))), stringsAsFactors=F, row.names=NULL)

cl_df <- data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",  "Type", "iter"))), stringsAsFactors=F,row.names=NULL)

cl_RFs <- list()

#For exporting models' results 
RF_results <- matrix(data=NA, nrow = 2, ncol = 5)

colnames(RF_results) <- c("Type", "Sens", "Spec", "AUC", "TSS")


#CV based on spatial unit
for (i in c(1:15)) {
  
  set.seed(i*754)
  ff <- unique(d_cutleish$Fold)[[i]]
  
  
  #data partitioning
  tra <- d_cutleish %>% filter(Fold != ff)
  tes <- d_cutleish %>% filter(Fold == ff)
  tra_mod <- as.data.frame(tra[names(tra) %in% n])
  tes_mod <- as.data.frame(tes[names(tes) %in% n])
  
  #tune
  o <- tune(cl_formula, data = tra_mod, rfq = TRUE)
  
  #model training
  cl_rf0_imb_rfq=imbalanced(cl_formula, 
                            ntree=500, 
                            data=tra_mod,
                            mtry = as.numeric(o$optimal[2]),
                            nodesize = as.numeric(o$optimal[1]),
                            method = "rfq",
                            do.trace=FALSE, 
                            importance="permute", 
                            statistics = T, 
                            imp.perf=T)
  cl_RFs <- c(cl_RFs, cl_rf0_imb_rfq)
  
  ##########
  
  #model validation and evaluation 
  test_out_0 <- predict(cl_rf0_imb_rfq, newdata=tes_mod)
  auc_out <- pROC::roc(response = test_out_0$yvar, predictor= test_out_0$predicted[,2], levels=c(0,1), auc = TRUE)
  best_threshold_out <- pROC::coords(auc_out, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_out <- best_threshold_out$sensitivity
  specificity_out <- best_threshold_out$specificity
  tss <- ecospat::ecospat.max.tss(test_out_0$predicted[,2], 
                                  tes_mod$Output)
  
  out_error <- data.frame(Type = 'cutaneous_leish', 
                          sens = max(sensitivity_out),
                          spec = max(specificity_out),
                          auc = auc_out$auc,
                          TSS = round(tss$max.TSS,2),
                          spatial_fold = i)
  oob_out_cv <- rbind(oob_out_cv, out_error)
  
  print(i)
}

cl_RFs_BM <- cl_RFs
cl_df_BM <- cl_df

### Performance
sens_performance <- aggregate(sens ~ Type, data = oob_out_cv,
                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                  lowerCI = Rmisc::CI(x, 0.95)[3],
                                                  upperCI = Rmisc::CI(x, 0.95)[1]))

spec_performance <- aggregate(spec ~ Type, data = oob_out_cv,
                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                  lowerCI = Rmisc::CI(x, 0.95)[3],
                                                  upperCI = Rmisc::CI(x, 0.95)[1]))

auc_performance <- aggregate(auc ~ Type, data = oob_out_cv,
                             FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                 lowerCI = Rmisc::CI(x, 0.95)[3],
                                                 upperCI = Rmisc::CI(x, 0.95)[1]))

tss_performance <- aggregate(TSS ~ Type, data = oob_out_cv,
                             FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                 lowerCI = Rmisc::CI(x, 0.95)[3],
                                                 upperCI = Rmisc::CI(x, 0.95)[1]))
#renaming performance measures 
sens_performance_BM <- sens_performance
spec_performance_BM <- spec_performance 
auc_performance_BM <- auc_performance
tss_performance_BM <- tss_performance

#saving model results in a final matrix
RF_results[1,1] <- "Base Model"
RF_results[1,2] <- round(sens_performance_BM$sens[1],2)
RF_results[1,3] <- round(spec_performance_BM$spec[1],2)
RF_results[1,4] <- round(auc_performance_BM$auc[1],2)
RF_results[1,5] <- round(tss_performance_BM$TSS[1],2)
RF_results

# end BASE MODEL


#### MACRO ALL MODEL ####

#creating formula
names(d_cutleish)
n <- setdiff(names(d_cutleish), c("Fold"))
cl_formula <- as.formula(paste("Output~", paste(n[-length(n)], collapse = " + ")))
cl_formula


##------------------------------------------------##
##   for loop for cutaneous leishmaniasis model   ##
##------------------------------------------------##

oob_out_cv <- data.frame(matrix(vector(),0, 6, dimnames=list(c(), c("auc", "sens", "spec", "TSS","Type", "spatial_fold"))),  stringsAsFactors=F, row.names=NULL)
cl_df <- data.frame(matrix(vector(), 0, 7, dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value", "Type", "iter"))), stringsAsFactors=F,row.names=NULL)

cl_RFs <- list()
importance_var <- list()

#CV based on spatial unit
for (i in c(1:15)) {
  
  set.seed(i*754)
  ff <- unique(d_cutleish$Fold)[[i]]
  
  ###########
  #data partitioning
  tra <- d_cutleish %>% filter(Fold != ff)
  tes <- d_cutleish %>% filter(Fold == ff)
  tra_mod <- as.data.frame(tra[names(tra) %in% n])
  tes_mod <- as.data.frame(tes[names(tes) %in% n])
  
  #tune
  o <- tune(cl_formula, data = tra_mod, rfq = TRUE)
  
  #model training
  cl_rf0_imb_rfq=imbalanced(cl_formula, 
                            ntree=500, 
                            data=tra_mod,
                            mtry = as.numeric(o$optimal[2]),
                            nodesize = as.numeric(o$optimal[1]),
                            method = "rfq",
                            do.trace=FALSE, 
                            importance="permute", 
                            statistics = T, 
                            imp.perf=T)
  cl_RFs <- c(cl_RFs, cl_rf0_imb_rfq)
  
  
  ##########
  #for PDP
  var <- colnames(tra_mod[-(ncol(tra_mod))])
  cl_df_1 <- data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), c("prob_1", "prob_se", "prob_adj", "variable", "variable_value",
                                                    "Type", "iter"))), stringsAsFactors=F,
                        row.names=NULL)
  
  for(j in 1:length(var)){
    hfp <- plot.variable(cl_rf0_imb_rfq, xvar.names = var[j], partial = TRUE,target = "1",  
                         npts = 25, grid.resoluion = 1)
    cl_df_0 <- data.frame(prob_1 = hfp$pData[[c(1, 2)]],
                          prob_se = hfp$pData[[c(1, 3)]], 
                          prob_adj = (hfp$pData[[c(1, 2)]] - min(hfp$pData[[c(1, 2)]]))/(max(hfp$pData[[c(1, 2)]])-min(hfp$pData[[c(1, 2)]])),
                          variables = var[j],
                          variable_value = hfp$pData[[c(1, 5)]],
                          Type =  "cutaneous_leish",
                          iter = i)
    cl_df_1 <- rbind(cl_df_1, cl_df_0)
    cl_df <- rbind(cl_df, cl_df_1)
    importance_var[[i]] <- data.frame(Variable = rownames(cl_rf0_imb_rfq$importance),
                                      Importance = abs(cl_rf0_imb_rfq$importance[,1]))
    
  }
  #end for PDP
  ############
  
  #model validation and evaluation
  test_out_0 <- predict(cl_rf0_imb_rfq, newdata=tes_mod)
  auc_out <- pROC::roc(response = test_out_0$yvar, predictor= test_out_0$predicted[,1], levels=c(0, 1), auc = TRUE)
  best_threshold_out <- pROC::coords(auc_out, "best", ret = c("threshold", "sensitivity", "specificity"))
  sensitivity_out <- best_threshold_out$sensitivity
  specificity_out <- best_threshold_out$specificity
  tss <- ecospat::ecospat.max.tss(test_out_0$predicted[,2], 
                                  tes_mod$Output)
  
  out_error <- data.frame(Type = 'cutaneous_leish', 
                          
                          sens = max(sensitivity_out),
                          spec = max(specificity_out),
                          auc = auc_out$auc,
                          TSS = round(tss$max.TSS,2),
                          spatial_fold = i)
  oob_out_cv <- rbind(oob_out_cv, out_error)
  
  print(i)
  
  
}

cl_RFs_MacroALL <- cl_RFs
cl_df_MacroALL <- cl_df

### Performance
sens_performance <- aggregate(sens ~ Type, data = oob_out_cv,
                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                  lowerCI = Rmisc::CI(x, 0.95)[3],
                                                  upperCI = Rmisc::CI(x, 0.95)[1]))


spec_performance <- aggregate(spec ~ Type, data = oob_out_cv,
                              FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                  lowerCI = Rmisc::CI(x, 0.95)[3],
                                                  upperCI = Rmisc::CI(x, 0.95)[1]))

auc_performance <- aggregate(auc ~ Type, data = oob_out_cv,
                             FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                 lowerCI = Rmisc::CI(x, 0.95)[3],
                                                 upperCI = Rmisc::CI(x, 0.95)[1]))

tss_performance <- aggregate(TSS ~ Type, data = oob_out_cv,
                             FUN = function(x) c(mean = Rmisc::CI(x, 0.95)[2],
                                                 lowerCI = Rmisc::CI(x, 0.95)[3],
                                                 upperCI = Rmisc::CI(x, 0.95)[1]))

sens_performance_MacroALL <- sens_performance
spec_performance_MacroALL <- spec_performance 
auc_performance_MacroALL <- auc_performance
tss_performance_MacroALL <- tss_performance

#result export
RF_results[2,1] <- "MACRO ALL"
RF_results[2,2] <- round(sens_performance$sens[1],2)
RF_results[2,3] <- round(spec_performance$spec[1],2)
RF_results[2,4] <- round(auc_performance$auc[1],2)
RF_results[2,5] <- round(tss_performance_MacroALL$TSS[1],2)
RF_results



################
### PD Plots ###
unique(cl_df_MacroALL$variables)
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Human_footprint"] <- "Human footprint"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Population"] <- "log(Population)"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Forest"] <- "Forest (cover)"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Annual_tmp"] <- "Temperature"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Annual_prec"] <- "Precipitation"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Annual_wet"] <- "Wet days"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Host_rich"] <- "Host richness"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Forest_los"] <- "Forest loss"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="Fun_rich"] <- "Functional richness"
cl_df_MacroALL$variables[cl_df_MacroALL$variables=="EVI_varian"] <- "EVI"


list_variables <- unique(cl_df_MacroALL$variables)

plot_list <- list()
for (i in 1:length(list_variables)) {
  if (i > 8) {
    vv <- list_variables[[i]]
    
    p <- ggplot(subset(cl_df_MacroALL, variables == vv),
                aes(x = variable_value, y = prob_1)) + 
      stat_smooth(aes(group = iter), 
                  color = 'lightgrey', method = 'loess', linewidth = 0.5, se = FALSE) + 
      stat_smooth(aes(), method = 'loess', linewidth = 1, se = FALSE, color = "#9ab522") +
      ylab(' ') +
      xlab(vv) +
      theme_bw(base_size = 12) + 
      theme(
        legend.position = "none",
        panel.grid = element_blank()  
      )
    plot_list[[i]] <- p
    
  } else {
    vv <- list_variables[[i]]
    
    p <- ggplot(subset(cl_df_MacroALL, variables == vv),
                aes(x = variable_value, y = prob_1)) + 
      stat_smooth(aes(group = iter), 
                  color = 'lightgrey', method = 'loess', linewidth = 0.5, se = FALSE) + 
      stat_smooth(aes(), method = 'loess', linewidth = 1, se = FALSE, color = "#14213d") +
      ylab(' ') +
      xlab(vv) +
      theme_bw(base_size = 12) + 
      theme(
        legend.position = "none",
        panel.grid = element_blank()  
      )
    plot_list[[i]] <- p
  }
}

pdp_all_thr <- grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)

text_grob <- textGrob("Probability",
                      x = unit(0.4, "npc"), 
                      y = unit(0.55, "npc"), 
                      rot = 90, 
                      gp = gpar(fontface = "bold", fontsize = 12))
pdp_all <- arrangeGrob(pdp_all_thr, left = text_grob)
grid.newpage() 
grid.draw(pdp_all)

#ggsave(file="PDP_MacroALL.png", pdp_all, width = 9, height = 6, units = "in", dpi=300)


###############
### Plot relative importance
importance_df <- bind_rows(importance_var)

# averaging
average_importance_df <- importance_df %>% group_by(Variable) %>% summarise(AverageImportance = mean(Importance))

# rel impo %
average_importance_df$relImp <- (average_importance_df$AverageImportance/
                                   (sum(average_importance_df$AverageImportance))*100)

# rearrange 
average_importance_df <- average_importance_df %>%
  arrange(desc(relImp))

####
average_importance_df$Variable[average_importance_df$Variable=="Human_footprint"] <- "Human footprint"
average_importance_df$Variable[average_importance_df$Variable=="Population"] <- "log(Population)"
average_importance_df$Variable[average_importance_df$Variable=="Forest"] <- "Forest (cover)"
average_importance_df$Variable[average_importance_df$Variable=="Annual_tmp"] <- "Temperature"
average_importance_df$Variable[average_importance_df$Variable=="Annual_pre"] <- "Precipitation"
average_importance_df$Variable[average_importance_df$Variable=="Annual_wet"] <- "Wet days"
average_importance_df$Variable[average_importance_df$Variable=="Host_rich"] <- "Host richness"
average_importance_df$Variable[average_importance_df$Variable=="EVI_varian"] <- "Enhanced Vegetation Index"
average_importance_df$Variable[average_importance_df$Variable=="BII"] <- "Biodiversity Intactness Index"
average_importance_df$Variable[average_importance_df$Variable=="Forest_los"] <- "Forest loss"
average_importance_df$Variable[average_importance_df$Variable=="Fun_rich"] <- "Functional richness"


average_importance_df$color <- ifelse(average_importance_df$Variable %in% 
                                        c("Biodiversity Intactness Index", 
                                          "Enhanced Vegetation Index", 
                                          "Functional richness",
                                          "Forest loss", 
                                          "Host richness"
                                        ), "#9ab522","#14213d")
#barplot

ggplot(average_importance_df, aes(x = reorder(Variable, relImp), y = relImp, fill = color)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  coord_flip() + 
  labs(y = "Relative importance (%)", x = "") +
  scale_fill_manual(values = c("#14213d", "#9ab522" )) +
  theme_bw(base_size = 14) + 
  theme(axis.ticks = element_line(color = "black"), 
        legend.position = "none")


#end MACRO ALL
###############



