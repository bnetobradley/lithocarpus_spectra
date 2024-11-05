## code for plotting plantnet classification of lithocarpus herbarium specimens 
## last updated by B Neto-Bradley: October 2024

## read in libraries
library(tidyverse)

## read in plantnet classification results
pn <- read_csv("data/plantnet_results_test1.csv")

## summarize plant net classification results
length(which(pn$Predicted == pn$Truth))/102
conf2 <- table(list(predicted=pn$Predicted, observed=pn$Truth))
resultsc <- caret::confusionMatrix(conf2)
results <- yardstick::conf_mat(conf2)
autoplot(results, type = "heatmap")
taxa <- unique(pn$Truth)

## plot pretty results
plantnet_cm <- ggplot(as.data.frame(resultsc$table), aes(x = predicted, y = observed, fill = Freq/6)) + geom_tile(color = "white") + geom_text(aes(label = round(Freq/6, digits = 2)), color = "black", size = 2) + coord_fixed() + scale_fill_gradient(low = "#e8dde8",high = "#895C87")  + scale_y_discrete(limits = rev(taxa)) + scale_x_discrete(labels=function(x) gsub("_", " ", x, fixed=TRUE), limits = taxa, position="top")+  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0), axis.text.y = element_text(hjust = 1, vjust = 0) ) + labs(x = "Pl@ntNet prediction from herbarium specimen images", y = "Ground Truth") + labs(fill = "") + theme(axis.text = element_text(face = 'italic'))
ggsave(filename = "figures/x_plantnet_confusionmatrix.png", plot = plantnet_cm, width = 15, height = 15,dpi = 600, units = "cm")


#### SPECTRA VS IMAGE BASED CLASSIFICATION
#### Author: Barbara Neto-Bradley
#### Last updated: October 2024

## the code below takes the spectral data previously cleaned in "1_cleaning_spectra.R"; for a set number of spectral scans (going to use 30, based on the work in 2.1), and looks at how well we can classify taxa, the same is done based on herbarium specimen images

## READ IN LIBRARIES
library(caret)
library(matrixStats)
library(pls)
library(spectrolab)
library(tidyverse)
library(MASS)
library(jsonlite)

## READ IN DATA
meta.df <- read_rds("data/clean_lithocapurs_metadat.rds")
spectra <- read_rds("data/clean_normalized_lithocarpus_spectra.rds")
spectra <- scale(spectra)

## REMOVE BLANK HERBARIUM PAPER & CLIP SCANS
meta.df <- meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip")
## remove multiple scans from same herbarium specimen, for now
meta.df <- meta.df[!duplicated(meta.df$Barcode),]

meta.df$species <- gsub("dachystachyus","dasystachyus",meta.df$species)
meta.df$species <- gsub("niewenhuisii","nieuwenhuisii",meta.df$species)

## summary data
summary <-  meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip") %>% group_by(species)  %>% summarise(n_distinct(samplename))

## KEEP ONLY BEST REPRESENTED SPECIES
sp <- c("leptogyne","gracilis", "nieuwenhuisii", "conocarpus", "coopertus", "ewyckii", "caudatifolius","urceolaris","dasystachyus","bennettii","sundaicus","elegans","luteus", "cantleyanus","lucidus","havilandii","echinifer")

meta.df <- meta.df %>% filter(species %in% sp)

## set seed for repeatability
set.seed(44)

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

## generate all unique combinations
c17 <- as.data.frame(t(combn(sp, 17)))
results <- c() 
spp <- c(c17[1,])

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


## storing specimen barcodes for images
test.barcodes<- test.meta.df[training.samples,]
train.barcodes<- test.meta.df[-training.samples,]
#write_csv(test.barcodes[c(1:20)],"test_barcodes_plantnet.csv")
#write_csv(train.barcodes[c(1:20)],"train_barcodes_plantnet.csv")

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
tocomp <- as.data.frame(cbind(train.barcodes$Barcode,test.transformed$species))
tocomp$predicted <- as.data.frame(predictions$class)
conf2 <- table(list(predicted=predictions$class, observed=test.transformed$species))
results <- caret::confusionMatrix(conf2)
# class.results <- as.data.frame(t(as.matrix(results$byClass)))
#class.results$obs <- paste(spp)
#results <- as.data.frame(t(as.matrix(results$overall)))
#results$obs <- paste(spp)

# big_results <- rbind(big_results, results)
#big_class <- rbind(big_class, class.results)

# b <- Sys.time()
#print(a-b)
#print(length(c17$V1)-i)


colnames(tocomp) <- c("s_Image","s_Truth","s_Predicted") 

## read in file with pn data
pn <- read_csv("data/plantnet_results_test1.csv")

pn <- pn %>% dplyr::select("Image","Truth","Predicted")
colnames(pn) <- c("p_Image","p_Truth","p_Predicted") 

pn$p_Image <- gsub(".jpg","",pn$p_Image)
pn$p_Truth <- gsub(" ","",pn$p_Truth)
pn$p_Predicted <- gsub(" ","",pn$p_Predicted)

bothmods <- as.data.frame(inner_join(pn, tocomp, by = c("p_Image"="s_Image")))
bothmods$disagree <- bothmods$p_Predicted != bothmods$s_Predicted

## calculate likelihoods of accurate ID GIVEN SPECIES

# Create an empty dataframe to store accuracy rates
likelihood_df <- data.frame(
  Category = character(),
  p_Accuracy = numeric(),
  s_Accuracy = numeric(),
  stringsAsFactors = FALSE
)

# Loop through unique values of p_Truth
for (i in 1:length(unique(bothmods$p_Truth))) {
  ii <- unique(bothmods$p_Truth)[i]
  
  # For p_Predicted and p_Truth
  p_total <- nrow(bothmods %>% filter(p_Predicted == ii))
  if (p_total > 0) {
    p_correct <- nrow(bothmods %>% filter(p_Predicted == ii & p_Truth == ii))
    p_accuracy <- p_correct / p_total  # Calculate accuracy
  } else {
    p_accuracy <- NA  # If no rows for the prediction, set to NA
  }
  
  # For s_Predicted and s_Truth
  s_total <- nrow(bothmods %>% filter(s_Predicted == ii))
  if (s_total > 0) {
    s_correct <- nrow(bothmods %>% filter(s_Predicted == ii & s_Truth == ii))
    s_accuracy <- s_correct / s_total  # Calculate accuracy
  } else {
    s_accuracy <- NA  # If no rows for the prediction, set to NA
  }
  
  # Add the results to the dataframe
  likelihood_df <- rbind(likelihood_df, data.frame(
    Category = ii,
    p_Accuracy = p_accuracy,
    s_Accuracy = s_accuracy
  ))
}

#write_csv(likelihood_df, "data/likelihood_df.csv")

### Plot

taxa <- unique(metspec$species)


lda_cm <- ggplot(as.data.frame(results$table), aes(x = predicted, y = observed, fill = Freq/6)) + geom_tile(color = "white") + geom_text(aes(label = round(Freq/6, digits = 2)), color = "black", size = 2) + coord_fixed() + scale_fill_gradient(low = "#e8dde8",high = "#895C87")  + scale_y_discrete(limits = rev(taxa)) + scale_x_discrete(labels=function(x) gsub("_", " ", x, fixed=TRUE), limits = taxa, position="top")+  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0), axis.text.y = element_text(hjust = 1, vjust = 0) ) + labs(x = "LDA prediction from herbarium specimen spectra", y = "Ground Truth") + labs(fill = "") + theme(axis.text = element_text(face = 'italic'))
ggsave(filename = "figures/x_lda_confusionmatrix.png", plot = lda_cm, width = 15, height = 15,dpi = 600, units = "cm")


## predict second test dataset

here_temp <- read.csv("data/joint_pn_spectra_validation.csv")
meta.df <- read_rds("data/clean_lithocapurs_metadat.rds")
spectra <- read_rds("data/clean_normalized_lithocarpus_spectra.rds")
spectra <- scale(spectra)

## REMOVE BLANK HERBARIUM PAPER & CLIP SCANS
meta.df <- meta.df %>% filter(unique_leaf != "page" & unique_leaf != "nothing" & unique_leaf != "??") %>% filter(! is.na(unique_leaf)) %>% filter(page_clip != "clip")
## remove multiple scans from same herbarium specimen, for now
meta.df <- meta.df[!duplicated(meta.df$Barcode),]
meta.df$species <- gsub("dachystachyus","dasystachyus",meta.df$species)
meta.df$species <- gsub("niewenhuisii","nieuwenhuisii",meta.df$species)
spec.df <- as.data.frame(as.vector(rownames(spectra)))
spec.df$hyperspectra <- as.matrix(spectra)
colnames(spec.df) <- c("sample_name", "hyperspectra")
metspec <- inner_join(meta.df, spec.df, by = c("samplename" = "sample_name"))

##
barcodes <- here_temp$Barcode
metspec <- metspec[metspec$Barcode %in% barcodes,]
test.data <- metspec[c(11,21)]

# Transform the data using the parameters from before
test.transformed <- preproc.param %>% predict(test.data)

# Make predictions
predictions <- predict(model, test.transformed$hyperspectra)
# Model accuracy
mean(predictions$class==test.transformed$species)
tocomp <- as.data.frame(cbind(metspec$Barcode,test.transformed$species))
tocomp$predicted <- as.data.frame(predictions$class)
#write.csv(tocomp,"data/spectra_joint_approach_results.csv")
conf2 <- table(list(predicted=predictions$class, observed=test.transformed$species))
results <- caret::confusionMatrix(conf2)

pnjson <- read_json("~/Downloads/all_predictions.json", simplifyVector = T)
workhere <- as.data.frame(pnjson$label_probabilities)
max_column <- apply(workhere, 1, function(row) {
  return(colnames(workhere)[which.max(row)])
})

plantnet_validation <- as.data.frame(cbind(pnjson$image_name,max_column))
plantnet_validation$V1 <- gsub(".jpg","",plantnet_validation$V1)
colnames(plantnet_validation)<-c("barcode","plantnet_prediction")
tocomp <- inner_join(plantnet_validation, tocomp, by = c("barcode"="V1"))
colnames(tocomp) <- c("barcode","plantnet_prediction","truth","spectra_predicted")

## algo to get the likelihood of accurate ID for each, select greater, return species name
spred <- as.data.frame(tocomp$spectra_predicted)

for (i in 2:102) {
  thisone <- which(likelihood_df$Category == tocomp$plantnet_prediction[i])
  plik <- likelihood_df$p_Accuracy[thisone]  
  thisone2 <- which(likelihood_df$Category == spred$`predictions$class`[i])
  slik <- likelihood_df$s_Accuracy[thisone2]
  lik <- cbind(plik,slik)
  all_lik <- rbind(all_lik,lik)
}

tocomp <- cbind(tocomp, as.data.frame(all_lik))
max <- apply(tocomp, 1, function(row){
  return(colnames(tocomp)[which.max(row)])
})
tocomp$max <- max

test <- c()
for (i in 1:102) {
  if(tocomp$max[i] == "slik"){
    test[i] <- print(as.character(tocomp$spectra_predicted$`predictions$class`[i]))
  } else {
    test[i] <- print(tocomp$plantnet_prediction[i]) }
}

tocomp$joint <- test

jointagreement <- c()
tocomp$jointagreement <- tocomp$plantnet_prediction == tocomp$spectra_predicted  
for (i in 1:102) {
  if(tocomp$jointagreement[i] == "TRUE"){
    jointagreement[i] <- print(tocomp$plantnet_prediction[i])
  } else {
    jointagreement[i] <- print(tocomp$joint[i])
  }
}

tocomp$jointagreement <- jointagreement

