#### PHYLOGENIES AND SPECTRA
#### Author: Barbara Neto-Bradley
#### Last updated: August 20th 2024
## the code below takes the spectral data previously cleaned in "1_cleaning_spectra.R"; for a set number of spectral scans (going to use 30-40, based on the work in 2.1), to ask if closely related taxa are more frequently classified based on their spectra (i.e. are closely related taxa also more closely related in spectral space)

## READ IN LIBRARIES
library(caret)
library(matrixStats)
library(pls)
library(spectrolab)
library(tidyverse)
library(ape)
library(phytools)
library(MASS)
library(yardstick)
## READ IN DATA
meta.df <- read_rds("data/clean_lithocapurs_metadat.rds")
spectra <- read_rds("data/clean_normalized_lithocarpus_spectra.rds")
spectra <- scale(spectra)
##phy <- read.tree("~/Documents/phd_thesis/2_phys_phy_fagaceae/data/input/supertimetree.hyboakconstrained.climateTips.tre")
phy <-read.tree("~/Downloads/Re_ 2018 Lithocarpus phylogeny/Faga_atpbrbcl2013_3sptree.phy")


## REMOVE BLANK HERBARIUM PAPER & CLIP SCANS
meta.df <- meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf !="??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip")
## remove multiple scans from same herbarium specimen, for now
meta.df <- meta.df[!duplicated(meta.df$Barcode),]
## update names for nieuwenhuisii & dasystachus
meta.df$species <- gsub("niewenhuisii", "nieuwenhuisii",meta.df$species) 
meta.df$species <-gsub("dachystachyus", "dasystachyus",meta.df$species) 

## FILTER SPECTRA to match records in meta.df 
spec.df <- as.data.frame(as.vector(rownames(spectra)))
spec.df$hyperspectra <- as.matrix(spectra)
colnames(spec.df) <- c("sample_name", "hyperspectra")
spd <- inner_join(spec.df, meta.df, by = c("sample_name" = "samplename"))
spd$binomial <- paste0("Lithocarpus_",spd$species)
phy$tip.label <- gsub("LITHOCARPUS", "Lithocarpus_",phy$tip.label)

## KEEP TAXA ON THE PHYLOGENY that have at least 30 spectral scans
spd <- spd[which(spd$binomial %in% phy$tip.label),]
keep_these <- spd %>% group_by(binomial) %>% summarise(which(n() > 30))
spd <- spd[which(spd$binomial %in% keep_these$binomial),]
spd <- spd %>% group_by(binomial) %>% sample_n(30)
## PRUNE THE TREE to keep only taxa with spectral data
phy <- keep.tip(phy, keep_these$binomial)
plot(phy)

### NEW CODE WITH LDA APPROACH
spec_mean_df <- spd[c(2,22)]
repeat{
training.samples <- spec_mean_df$binomial %>%
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
model <- lda(x = train.transformed$hyperspectra, grouping = train.transformed$binomial)
# Make predictions
predictions <- predict(model, test.transformed$hyperspectra)
addthis <- as.data.frame(predictions$class)
big_predictions <- rbind(big_predictions, addthis)
if(nrow(big_predictions)>7800){break}}

## add in truth
big_predictions$truth <- test.transformed$binomial

# Model accuracy
mean(big_predictions$`predictions$class`== big_predictions$truth)

conf2 <- table(list(predicted=big_predictions$`predictions$class`, observed=big_predictions$truth))
resultsc <- caret::confusionMatrix(conf2)
results <- conf_mat(conf2)
autoplot(results, type = "heatmap")

## for each pair, at what rate is spA misclassified as sp B
## error rates between species pairs
tudo <- as.data.frame(results$table)

## get MRCA between each taxon pair
tudo$node <- NA
ac_datalist = list()
for (i in 1:length(tudo$predicted)) {
  this_pair <- c(paste(tudo[i,1]), paste(tudo[i,2]))
  ac_datalist[i] <- getMRCA(phy = phy, tip = this_pair)
  tudo$node[i] <- ac_datalist[i] }
tudo$node <- as.numeric(tudo$node)

## bind time since MRCA to misclassification error rates
ac_ages <- as.data.frame(branching.times(phy))
colnames(ac_ages) <- "time_mya"
ac_ages$node <- rownames(ac_ages)
ac_ages$node <- as.numeric(ac_ages$node)
tudo <- full_join(tudo, ac_ages, by = c("node" = "node"))
tudo_plot <- tudo
tudo <- tudo[-which(tudo$predicted == tudo$observed),]
## plot error as a function of divergence time
mod <- lm(tudo$Freq ~ tudo$time_mya)
#tudo <- tudo[68,]
summary(mod)
anova(mod)
#png("figures/species_misclassification_mrca.png", units = "cm", width = 16, height = 10, res = 600)
 gb <- ggplot(tudo,aes(time_mya, Freq/606)) + geom_jitter(width = 0.25, height = 0.005, alpha = 0.5,color = "#895C87")+  theme_minimal() + labs(x = "Time since species' MRCA (MY)", y = "Times species A is misclassified as B") 
dev.off()

tudo$grp <- as.factor(tudo$time_mya)
ga <- ggplot(as.data.frame(resultsc$table), aes(x = predicted, y = observed, fill = Freq/606)) + geom_tile(color = "white") + geom_text(aes(label = round(Freq/606, digits = 2)), color = "black", size = 2) + coord_fixed() + scale_fill_gradient(low = "#e8dde8",high = "#895C87")  + scale_y_discrete(limits = phy$tip.label) + scale_x_discrete(labels=function(x) gsub("Lithocarpus_", " ", x, fixed=TRUE), limits = rev(phy$tip.label), position="top")+  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0), axis.text = element_text(face = 'italic'), axis.text.y = element_text(hjust = 1, vjust = 0) ) + labs(x = "Predicted", y = "Observed") + labs(fill = "")
  

plot.phylo(phy, use.edge.length = F)


ggsave(filename = "figures/3_heatmap_0824.png", plot = ga, width = 15, height = 15,dpi = 600, units = "cm")
ggsave(filename = "figures/3_mrca_misclassification_0824.png", plot = gb, width = 10, height = 10, dpi = 600, units = "cm")
