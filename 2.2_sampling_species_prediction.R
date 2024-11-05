#### HOW MANY SPECIES CAN WE DISTINGUISH BETWEEN IN ONE GO
#### Author: Barbara Neto-Bradley
#### Last updated: June 30th 2023

## the code below takes the spectral data previously cleaned in "1_cleaning_spectra.R"; for a set number of spectral scans (going to use 30-40, based on the work in 2.1), how does increasing the number of species present in the set affect our power to perceive them?


## need to generate list of species for which we have at least 30 specimens
## ideally would create a suite of unique combinations of these species
## but I think this will be computationally a lot. an alternative is picking x number of species at random for each 'bin' of choice and repeating this many times. you won't be able to discriminate the effect of a particular single species, but it will give an overall general sense of the change in discriminatory power
## will start with one order for these species, picked at random
## for each suite of unique combinations, need to test the accuracy of discriminating 2, 3, 4 ,..., n(species(specimens > 30))

## READ IN LIBRARIES
library(caret)
library(matrixStats)
library(pls)
library(spectrolab)
library(tidyverse)
library(MASS)

## READ IN DATA
## READ IN DATA
meta.df <- read_rds("data/clean_lithocapurs_metadat.rds")
spectra <- read_rds("data/clean_normalized_lithocarpus_spectra.rds")
spectra <- scale(spectra)

## REMOVE BLANK HERBARIUM PAPER & CLIP SCANS
meta.df <- meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip")
## remove multiple scans from same herbarium specimen, for now
meta.df <- meta.df[!duplicated(meta.df$Barcode),]

## summary data
summary <-  meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip") %>% group_by(species)  %>% summarise(n_distinct(samplename))

## KEEP ONLY BEST REPRESENTED SPECIES
sp <- c("leptogyne","gracilis", "niewenhuisii", "conocarpus", "coopertus", "ewyckii", "caudatifolius","urceolaris","dachystachyus","bennettii","sundaicus","elegans","luteus", "cantleyanus","lucidus","havilandii","echinifer")

meta.df <- meta.df %>% filter(species %in% sp)

## KEEP ONLY 30 SCANS PER SPECIES
meta.df <- meta.df %>% group_by(species) %>% sample_n(30)

## get list of species
sp <- unique(meta.df$species)

## make a place for results
spec.df <- as.data.frame(as.vector(rownames(spectra)))
spec.df$hyperspectra <- as.matrix(spectra)
colnames(spec.df) <- c("sample_name", "hyperspectra")
metspec <- inner_join(meta.df, spec.df, by = c("samplename" = "sample_name"))

#saveRDS(spec.df, file = "data/lithocarpus_spectra.rds")

## make a place for results
big_results <- matrix(nrow = 10)
  
## generate all unique combinations
c17 <- as.data.frame(t(combn(sp, 17)))

## create a sequence with increasing amounts of herbarium records per species eg. c(10,20,30,40,50,60,70,80,90,100); and test predictive efficacy
 
## replace 100 in line 66 with length(c9$V1)

results <- c() 
for (i in 1:length(c17$V1)) {
  a <- Sys.time()
    spp <- c(c17[i,])
    ## KEEP given SPECIES
    ## KEEP i SCANS PER SPECIES
    test.meta.df <- metspec %>% filter(species %in% spp ) 
    
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
    model <- lda(x = train.transformed$hyperspectra, grouping = train.transformed$species)
    # Make predictions
    predictions <- predict(model, test.transformed$hyperspectra)
    # Model accuracy
    mean(predictions$class==test.transformed$species)
    
    conf2 <- table(list(predicted=predictions$class, observed=test.transformed$species))
    results <- caret::confusionMatrix(conf2)
   # class.results <- as.data.frame(t(as.matrix(results$byClass)))
    #class.results$obs <- paste(spp)
    results <- as.data.frame(t(as.matrix(results$overall)))
    #results$obs <- paste(spp)
    
    big_results <- rbind(big_results, results)
    #big_class <- rbind(big_class, class.results)
    
    b <- Sys.time()
    print(a-b)
    print(length(c17$V1)-i)
}




#c9 <- cbind(big_results, c9)

###c10 <- cbind(big_results, c10)
###c8 <- cbind(big_results, c8)
####c11 <- cbind(big_results, c11) 
####c7 <- cbind(big_results, c7)
####c12 <- cbind(big_results, c12) 
####c13 <- cbind(big_results, c13)
####c14 <- cbind(big_results, c14)
####c15 <- cbind(big_results, c15)
####c16 <- cbind(big_results, c16)
####c17 <- cbind(big_results, c17)
####c6 <- cbind(big_results, c6) 
####c5 <- cbind(big_results, c5)
####c4 <- cbind(big_results, c4)
####c3 <- cbind(big_results, c3)
####c2 <- cbind(big_results, c2)
##c17 <- cbind(results, c17)

saveRDS(c2, file = "data/results/0824_c2.rds")
saveRDS(c3, file = "data/results/0824_c3.rds")
saveRDS(c4, file = "data/results/0824_c4.rds")
saveRDS(c5, file = "data/results/0824_c5.rds")
saveRDS(c6, file = "data/results/0824_c6.rds")
saveRDS(c7, file = "data/results/0824_c7.rds")
saveRDS(c8, file = "data/results/0824_c8.rds")
saveRDS(c9, file = "data/results/0824_c9.rds")
saveRDS(c10, file = "data/results/0824_c10.rds")
saveRDS(c11, file = "data/results/0824_c11.rds")
saveRDS(c12, file = "data/results/0824_c12.rds")
saveRDS(c13, file = "data/results/0824_c13.rds")
saveRDS(c14, file = "data/results/0824_c14.rds")
saveRDS(c15, file = "data/results/0824_c15.rds")
saveRDS(c16, file = "data/results/0824_c16.rds")
saveRDS(c17, file = "data/results/0824_c17.rds")

#### START HERE IF PLOTTING ####
c2 <- read_rds("data/results/0824_c2.rds")
c3 <- read_rds("data/results/0824_c3.rds")
c4 <- read_rds("data/results/0824_c4.rds")
c5 <- read_rds("data/results/0824_c5.rds")
c6 <- read_rds("data/results/0824_c6.rds")
c7 <- read_rds("data/results/0824_c7.rds")
c8 <- read_rds("data/results/0824_c8.rds")
c9 <- read_rds("data/results/0824_c9.rds")
c10 <- read_rds("data/results/0824_c10.rds")
c11 <- read_rds("data/results/0824_c11.rds")
c12 <- read_rds("data/results/0824_c12.rds")
c13 <- read_rds("data/results/0824_c13.rds")
c14 <- read_rds("data/results/0824_c14.rds")
c15 <- read_rds("data/results/0824_c15.rds")
c16 <- read_rds("data/results/0824_c16.rds")
c17 <- read_rds("data/results/0824_c17.rds")

##summarise
big_mean <- as.data.frame(c(mean(c2$Accuracy) , mean(c3$Accuracy), mean(c4$Accuracy), mean(c5$Accuracy), mean(c6$Accuracy), mean(c7$Accuracy), mean(c8$Accuracy), mean(c9$Accuracy), mean(c10$Accuracy), mean(c11$Accuracy),  mean(c12$Accuracy), mean(c13$Accuracy), mean(c14$Accuracy), mean(c15$Accuracy), mean(c16$Accuracy),mean(c17$Accuracy)))
big_mean$obs <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
big_median <- as.data.frame(c(median(c2$Accuracy) , median(c3$Accuracy), median(c4$Accuracy), median(c5$Accuracy), median(c6$Accuracy),  median(c7$Accuracy),median(c8$Accuracy), median(c9$Accuracy),median(c10$Accuracy),median(c11$Accuracy),median(c12$Accuracy), median(c13$Accuracy), median(c14$Accuracy), median(c15$Accuracy), median(c16$Accuracy),median(c17$Accuracy)))
big_median$obs <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

big_sd <- as.data.frame(c(sd(c2$Accuracy) , sd(c3$Accuracy), sd(c4$Accuracy), sd(c5$Accuracy), sd(c6$Accuracy),sd(c7$Accuracy),sd(c8$Accuracy),sd(c9$Accuracy),sd(c10$Accuracy),sd(c11$Accuracy),  sd(c12$Accuracy), sd(c13$Accuracy), sd(c14$Accuracy), sd(c15$Accuracy), sd(c16$Accuracy),sd(c17$Accuracy)))
big_sd$obs <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
colnames(big_mean) <- c("V1","V2")
big <- cbind(big_mean,big_sd, big_median)
colnames(big) <- c("mean","species","sd","speciess", "median", "speciesss")

## plot
exportthisplot <- ggplot(NULL, aes(x = results)) +
  theme_classic() + xlab("Number of species") + ylab("Accuracy") +
  #geom_boxplot(data = big_result_s, position = "dodge2", width = 5, fill = NA, colour = "#284e7c", outlier.alpha = 0, alpha = 0.25) +
  geom_jitter(data = c2, aes(y = Accuracy, x = 2), width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c3, aes(y = Accuracy, x = 3),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c4, aes(y = Accuracy, x = 4),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c5, aes(y = Accuracy, x = 5),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c6, aes(y = Accuracy, x = 6),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c7, aes(y = Accuracy, x = 7),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c8, aes(y = Accuracy, x = 8),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c9, aes(y = Accuracy, x = 9),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c10, aes(y = Accuracy, x = 10),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15)+
  geom_jitter(data = c11, aes(y = Accuracy, x = 11),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15)+
  geom_jitter(data = c12, aes(y = Accuracy, x = 12), width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c13, aes(y = Accuracy, x = 13),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c14, aes(y = Accuracy, x = 14),width = 0.3, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c15, aes(y = Accuracy, x = 15), width = 0.3,size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c16, aes(y = Accuracy, x = 16), width = 0,size = 0.5, color = "#BDA3BC", alpha = 0.15)+
  geom_jitter(data = c17, aes(y = Accuracy, x = 17), width = 0,size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_point(data = big_mean, aes(y = V1, x = V2), shape = 8, colour = "#895C87") +
  geom_errorbar(data = big, aes(x = big$species, ymin = big$mean - big$sd, ymax = big$mean + big$sd), colour = "#895C87", width = 0.5) +
  scale_x_continuous(breaks = c(2:17), labels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1")) 

## save plot in set dimensions
ggsave(filename = "figures/2_accuracy_numberofspecies_0824.png", plot = exportthisplot, width = 12, height = 10, dpi = 600, units = "cm")

## putting together accuracy and kappa across all iterations

brs <- rbind(c2[,1:2],c3[,1:2],c4[,1:2],c5[,1:2],c6[,1:2],c7[,1:2],c8[,1:2],c9[,1:2],c10[,1:2],c11[,1:2],c12[,1:2],c13[,1:2],c14[,1:2],c15[,1:2],c16[,1:2],c17[,1:2])

brs$species <- c(rep.int(2,136),rep.int(3,680),rep.int(4,2380),rep.int(5,6188),rep.int(6,12376),rep.int(7,19448),rep.int(8,24310),rep.int(9,24310),rep.int(10,19448),rep.int(11,12376),rep.int(12,6188),rep.int(13,2380),rep.int(14,680),rep.int(15,136),rep.int(16,17),rep.int(17,1))

mod <- lm(data = brs, Accuracy~species)
anova(mod)
summary(mod)

moda <- lm(data = brs, Kappa~species)
anova(moda)
summary(moda)

bigk_mean <- as.data.frame(c(mean(c2$Kappa) , mean(c3$Kappa), mean(c4$Kappa), mean(c5$Kappa), mean(c6$Kappa), mean(c7$Kappa), mean(c8$Kappa), mean(c9$Kappa), mean(c10$Kappa), mean(c11$Kappa),  mean(c12$Kappa), mean(c13$Kappa), mean(c14$Kappa), mean(c15$Kappa), mean(c16$Kappa),mean(c17$Kappa)))

bigk_sd <- as.data.frame(c(sd(c2$Kappa) , sd(c3$Kappa), sd(c4$Kappa), sd(c5$Kappa), sd(c6$Kappa),sd(c7$Kappa),sd(c8$Kappa),sd(c9$Kappa),sd(c10$Kappa),sd(c11$Kappa),  sd(c12$Kappa), sd(c13$Kappa), sd(c14$Kappa), sd(c15$Kappa), sd(c16$Kappa),sd(c17$Kappa)))
bigk_sd$obs <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

bigk <- cbind(bigk_mean,bigk_sd)
colnames(bigk) <- c("mean","sd","obs")
## plot kappa in same style as accuracy
exportthiskplot <- ggplot(NULL, aes(x = results)) +
  theme_classic() + xlab("Number of species") + ylab("Kappa") +
  #geom_boxplot(data = big_result_s, position = "dodge2", width = 5, fill = NA, colour = "#284e7c", outlier.alpha = 0, alpha = 0.25) +
  geom_jitter(data = c2, aes(y = Kappa, x = 2), width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c3, aes(y = Kappa, x = 3),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c4, aes(y = Kappa, x = 4),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c5, aes(y = Kappa, x = 5),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c6, aes(y = Kappa, x = 6),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c7, aes(y = Kappa, x = 7),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c8, aes(y = Kappa, x = 8),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c9, aes(y = Kappa, x = 9),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c10, aes(y = Kappa, x = 10),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15)+
  geom_jitter(data = c11, aes(y = Kappa, x = 11),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15)+
  geom_jitter(data = c12, aes(y = Kappa, x = 12), width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c13, aes(y = Kappa, x = 13),width = 0.3,height = 0.025, size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c14, aes(y = Kappa, x = 14),width = 0.3, height = 0.025,size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c15, aes(y = Kappa, x = 15), width = 0.3,height = 0.025,size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_jitter(data = c16, aes(y = Kappa, x = 16), width = 0,height = 0.025,size = 0.5, color = "#BDA3BC", alpha = 0.15)+
  geom_jitter(data = c17, aes(y = Kappa, x = 17), width = 0,height = 0.025,size = 0.5, color = "#BDA3BC", alpha = 0.15) +
  geom_point(data = bigk, aes(y = mean, x = obs), shape = 8, colour = "#895C87") +
  geom_errorbar(data = bigk, aes(x = bigk$obs, ymin = bigk$mean - bigk$sd, ymax = bigk$mean + bigk$sd), colour = "#895C87", width = 0.5) +
  scale_x_continuous(breaks = c(2:17), labels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17")) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c("0","0.2","0.4","0.6","0.8","1")) 

ggsave(filename = "figures/2_kappa_numberofspecies_0824.png", plot = exportthiskplot, width = 12, height = 10, dpi = 600, units = "cm")
