#### HOW MANY OBSERVATIONS DO WE NEED TO PREDICT SPECIES FROM SPECTRA
#### Author: Barbara Neto-Bradley
#### Last updated: August 6th 2024

## the code below takes the spectral data previously cleaned in "1_cleaning_spectra.R" and uses the three best represented species, to identify how many observations of herbarium specimens are needed per taxa to build a model which accurately predicts species identity
## n.b. skip to "B. START HERE TO PLOT DATA" if revisiting for cosmetic touch ups

## READ IN LIBRARIES
library(spectrolab)
library(plotrix)
library(caret)
library(matrixStats)
library(pls)
library(spectrolab)
library(tidyverse)

## A. START HERE TO RUN ANALYSIS ####
## READ IN DATA
meta.df <- read_rds("data/clean_lithocapurs_metadat.rds")
spectra <- read_rds("data/clean_normalized_lithocarpus_spectra.rds")
spectra <- scale(spectra)

## summary data
summary <-  meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip") %>% group_by(species)  %>% summarise(n_distinct(samplename))

## REMOVE BLANK HERBARIUM PAPER & CLIP SCANS
meta.df <- meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip")

## KEEP ONLY BEST REPRESENTED 3 SPECIES
## LEPTOGYNE, GRACILIS AND NIEWENHUISII
meta.df <- meta.df %>% filter(species == "leptogyne" | species == "gracilis" | species == "niewenhuisii")

## KEEP ONLY ONE SCAN PER HERBARIUM SHEET
meta.df <- meta.df[!duplicated(meta.df$Barcode),]

## MATCH UP DATASETS 
spec.df <- as.data.frame(as.vector(rownames(spectra)))
spec.df$hyperspectra <- as.matrix(spectra)
colnames(spec.df) <- c("sample_name", "hyperspectra")
metspec <- inner_join(meta.df, spec.df, by = c("samplename" = "sample_name"))

## make a place for results
big_results <- matrix(ncol = 7)
big_class <- matrix(ncol = 11)

## create a sequence with increasing amounts of herbarium records per species eg. c(10,20,30,40,50,60,70,80,90,100); and test predictive efficacy
start.time <- Sys.time()
repeat {
  for (i in c(10,20,30,40,50,60,70,80,90,100)) {
    
    ## KEEP i SCANS PER SPECIES
    lepto_i <- metspec %>% filter(species == "leptogyne") %>% sample_n(i)
    graci_i <- metspec %>% filter(species == "gracilis") %>% sample_n(i)
    niewe_i <- metspec %>% filter(species == "niewenhuisii") %>% sample_n(i)
    test.meta.df <- rbind(lepto_i, graci_i, niewe_i)
    
    ##
    barcodes <- test.meta.df$Barcode
    spec_mean_df <- test.meta.df[c(11,21)]
    training.samples <- spec_mean_df$species %>%
      createDataPartition(p = 0.8, list = FALSE)
    train.data <- spec_mean_df[training.samples, ]
    test.data <- spec_mean_df[-training.samples, ]
    
    # Estimate preprocessing parameters
    preproc.param <- train.data %>% 
      preProcess(method = c("center", "scale"))
    # Transform the data using the estimated parameters
    train.transformed <- preproc.param %>% predict(train.data)
    test.transformed <- preproc.param %>% predict(test.data)
    
    # Fit the model
    model <- lda(x = train.transformed[,2], grouping = train.transformed[,1])
    # Make predictions
    predictions <- predict(model, test.transformed[,2])
    # Model accuracy
    mean(predictions$class==test.transformed[,1])
    
    conf2 <- table(list(predicted=predictions$class, observed=test.transformed$species))
    results <- caret::confusionMatrix(conf2)
    class.results <- as.data.frame(as.matrix(results$byClass))
    class.results$obs <- paste(i)
    results <- as.data.frame(t(as.matrix(results$overall)))
    results$obs <- paste(i)
  
  big_results <- rbind(big_results, results)
  big_class <- rbind(big_class, class.results)
  }
  
  if((dim(big_results)[1]) > 1000) {
    break 
  }
  }
  end.time <- Sys.time()
  end.time - start.time

## plotting out results
  ggplot() + geom_jitter(aes(as.numeric(big_results$obs), big_results$Accuracy)) + theme_minimal()  
  
  
  
#### EVERYTHING HERE IS PLS-DA FOCUSED
## TO DO: EVALUATE HOW MANY COMPONENTS TO USE
#plsda.fit <- plsda(metspec$hyperspectra, metspec$species, ncomp = 30)
#perf.plsda <- perf(plsda.fit, validation = "Mfold", folds = 5, progressBar = TRUE, auc = TRUE, nrepeat = 50)

#plot(perf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
#plot(perf.plsda)
## select number of components
#perf.plsda$choice.ncomp

## make a place for results

#big_results <- matrix(ncol = 7)
## create a sequence with increasing amounts of herbarium records per species eg. c(10,20,30,40,50,60,70,80,90,100); and test predictive efficacy

#repeat {
 # results <- c() 
#for (i in c(10,20,30,40,50,60,70,80,90,100)) {

## KEEP i SCANS PER SPECIES
#lepto_i <- metspec %>% filter(species == "leptogyne") %>% sample_n(i)
#graci_i <- metspec %>% filter(species == "gracilis") %>% sample_n(i)
#niewe_i <- metspec %>% filter(species == "niewenhuisii") %>% sample_n(i)
#test.meta.df <- rbind(lepto_i, graci_i, niewe_i)

#partition the data: 70% for training, 30% for testing
#training.samples <- createDataPartition(test.meta.df$species, p = 0.7, list = FALSE)
#train.data  <- test.meta.df[training.samples, ]
#test.data <- test.meta.df[-training.samples, ]

## For PLS-DA, train the model
#plsda.train <- plsda(train.data$hyperspectra, train.data$species, ncomp = 20)
# then predict
#test.predict <- predict(plsda.train, test.data$hyperspectra, dist = "max.dist")
# store prediction for the 30th component
#prediction <- test.predict$class$max.dist[,20] 
# calculate the error rate of the model
#mat10 <- confusionMatrix(as.factor(prediction),as.factor(test.data$species))
#acc <- t(mat10$overall)
#results <- rbind(results, acc)
#}
#  big_results <- rbind(big_results, results)
  
#  if((dim(big_results)[1]) > 1000) {
#    break 
#    }
#}

## calculate average for each line in big_results
#big_results <- as.data.frame(big_results)
#big_results$V1 <- c(NA, rep.int(c(10,20,30,40,50,60,70,80,90,100),100))
#rownames(big_results) <- big_results$V1
#big_results <- big_results %>% select(-V1)
#big_results %>% pivot_wider()
#br <- as.matrix(big_results)
#rowMedians(br)
#big_result_s <- gather(big_results,"source_col","n_specimens", 1:100)
## save file so don't have to rerun
#saveRDS(file = "data/nrecordspredictionaccuracy.rds", big_result_s)
#saveRDS(file = "data/nrecordspredictionaccuracy_0824.rds", big_results)

## B. START HERE TO PLOT DATA####
## unhash below to read in file & pick up from plotting
  
#saveRDS(file = "data/nrecordspredictionaccuracy_0824.rds", big_results)
# big_results <- read_rds("data/nrecordspredictionaccuracy_0824.rds")
brs <- big_results %>% group_by(obs) %>% summarise(mean(Accuracy), median(Accuracy),sd(Accuracy),mean(Kappa), median(Kappa),sd(Kappa))

library(sjPlot)
library(segmented)

## fit simple linear regression model
fit <- lm(Accuracy~obs, data = big_results)
# use davies test to find initial breakpoint estimate
davies.test(fit)

## fit piecewise regression model to original model
segmented_model <- segmented(fit, seg.Z = ~obs,  psi = 30)

summary(segmented_model)
breakpoint <- segmented_model$psi[,"Est."]
print(breakpoint)
slope(segmented_model)

confint.segmented(segmented_model)

segmented.fit$coefficients
## compare simple linear model to breakpoint model
anova(fit,segmented.fit)
AIC(fit,segmented.fit)

## plot kapppa with mean and standard deviation overlaid
exportplot <- ggplot(NULL, aes(x = obs, group = obs)) +
   theme_classic() + xlab("Specimens Included Per Species") + ylab("Kappa") +
  #geom_boxplot(data = big_result_s, position = "dodge2", width = 5, fill = NA, colour = "#284e7c", outlier.alpha = 0, alpha = 0.25) +
  geom_jitter(data = big_results, aes(y = Kappa), size = 1, width = 2, color = "#BDA3BC", alpha = 0.5)  +
 geom_errorbar(data = brs, aes(x = brs$obs, ymin = brs$`mean(Kappa)`-brs$`sd(Kappa)`, ymax = brs$`mean(Kappa)` + brs$`sd(Kappa)`), colour = "#895C87", width = 3) +
  geom_point(data = brs, aes(y = `mean(Kappa)`, x = brs$obs), shape = 8, colour = "#895C87") +
  scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100), labels = c("10","20","30","40","50","60","70","80","90","100")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1")) 

## save plot in set dimensions
ggsave(filename = "figures/1_kappa_obsperspecies_0824.png", plot = exportplot, width = 12, height = 10, dpi = 600, units = "cm")

## plot accuracy with mean and standard deviation overlaid
exportplotacc <- ggplot(NULL, aes(x = obs, group = obs)) +
  theme_classic() + xlab("Specimens Included Per Species") + ylab("Accuracy") +
  #geom_boxplot(data = big_result_s, position = "dodge2", width = 5, fill = NA, colour = "#284e7c", outlier.alpha = 0, alpha = 0.25) +
  geom_jitter(data = big_results, aes(y = Accuracy), size = 1, width = 2, color = "#BDA3BC", alpha = 0.5)  +
  geom_errorbar(data = brs, aes(x = brs$obs, ymin = brs$`mean(Accuracy)`-brs$`sd(Accuracy)`, ymax = brs$`mean(Accuracy)` + brs$`sd(Accuracy)`), colour = "#895C87", width = 3) +
  geom_point(data = brs, aes(y = `mean(Accuracy)`, x = brs$obs), shape = 8, colour = "#895C87") +
  scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100), labels = c("10","20","30","40","50","60","70","80","90","100")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1")) 

## save plot in set dimensions
ggsave(filename = "figures/1_accuracy_obsperspecies_0824.png", plot = exportplotacc, width = 12, height = 10, dpi = 600, units = "cm")



## plot with breakpoint line
pred_seg <- as.data.frame(predict(segmented.fit, interval = "confidence"))
pred_seg$obss <- c(10,rep.int(c(10,20,30,40,50,60,70,80,90,100),100))
big_results <- cbind(big_results[1:1001,1:8],pred_seg)

breakplot <- ggplot(data = big_results, aes(x = obs)) +
  geom_jitter(data = big_results, aes(y = Accuracy), size = 1, width = 2, color = "#BDA3BC", alpha = 0.5)+
#geom_ribbon(aes(ymin = lwr, ymax = upr), fill ="#895C87", alpha = 0.15) +  #geom_line(aes( y = fit), size = 0.75,color = "#895C87") + #geom_abline(aes(intercept = 0.450556, slope = 0.009750)) +
  #(x*0.0008489 + 0.8300794) and (x*0.0065 + 0.6337037)
  geom_function(fun = Vectorize(function(x){
    if(x>31.728)
      return(NA)
    else
      return(x*0.0075876 + 0.6547948)
  }), color = '#895C87', cex = 1)  +
  geom_function(fun = Vectorize(function(x){
    if(x<31.728)
      return(NA)
    else
      return(x*0.0004643 + 0.8782211)
  }), color = '#895C87', cex =1) +
#geom_abline(aes(intercept = 0.7451317, slope = 0.001273)) +
  theme_classic() + xlab("Specimens Included Per Species") + ylab("Accuracy") +  scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100), labels = c("10","20","30","40","50","60","70","80","90","100")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1")) 

## save plot in set dimensions
ggsave(filename = "figures/1_accuracy_breakplot_obsperspecies_0824.png", plot = breakplot, width = 12, height = 10, dpi = 600, units = "cm")
