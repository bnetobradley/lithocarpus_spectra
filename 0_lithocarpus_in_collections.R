## 28st October 2024
## Exploration of GBIF data

## load libraries
library(tidyverse)
## GBIF.org (28 October 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.7a4m94

## read in raw gbif data (in darwin core format
##*old* lith <- read.delim(file = "~/Documents/phd_thesis/placeholder_until_gitsetup/data_exploration/0033011-210914110416597/occurrence.txt")
lith <- read.delim(file = "data/0008178-241024112534372/occurrence.txt")

## pull out columns of interest: species name, date of record, collection, location of collection, basis of record

lithsum <- lith %>% select(collectionCode, basisOfRecord, reproductiveCondition, occurrenceStatus, year, month, day, countryCode,acceptedScientificName, publishingCountry, iucnRedListCategory, institutionCode)

nrecords <- lithsum %>% group_by(acceptedScientificName) %>% summarise(n())

## make a list of the number of records per species
nspecimens <- lithsum %>% filter(basisOfRecord == "PRESERVED_SPECIMEN") %>% group_by(acceptedScientificName) %>% summarise(n())

## make a list of species with iucn status
iucnspecimens <- lithsum %>% filter(iucnRedListCategory != "") %>% select(acceptedScientificName, iucnRedListCategory) %>% unique()

## combine the two lists into one
lithocarpus_list <- full_join(nspecimens, iucnspecimens)

xl <- readxl::read_xlsx("~/Documents/phd_thesis/placeholder_until_gitsetup/lithocarpus.xlsx")
xl$speciesmatch <- paste(xl$Genus, xl$Species, xl$Authorship)

new_lithocarpus <- full_join(lithocarpus_list, xl, by = c("acceptedScientificName" =  "speciesmatch"))
#write.csv(new_lithocarpus, file = "lith_list_records.csv")


## summarise how many records are listed only to genus
countries <- lithsum %>% filter(acceptedScientificName == "Lithocarpus Blume") %>% filter(basisOfRecord == "PRESERVED_SPECIMEN") %>% group_by(institutionCode) %>% summarise(n())

## temporal bias in unidentified records
unid <- lithsum %>% filter(acceptedScientificName != "Lithocarpus Blume")

## layer histograms of records proportion of un-ID'd
ggplot() + geom_histogram(data = lithsum, aes(x =year, fill = "Records identified to genus"), binwidth = 10, color = "#FFFFFF") + geom_histogram(data = unid, aes(x = year, fill = "Records identified to species"), binwidth = 10, color = "#FFFFFF") + theme_classic() + xlim(c(1850, 2021)) + xlab("Year of Collection") +  ylab("Number of Records Identified") + theme(axis.text=element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=12), legend.title = element_blank(), legend.position = c(0.3, 0.9)) + scale_fill_manual(values = c("#895C87","#B599B4"))


# Summarize counts of genus and species by year bin
lithsum_summary <- lithsum %>% mutate(year_bin = floor(year / 5) * 5) %>%  group_by(year_bin) %>% summarize(count_genus = n())  # Count genus records per bin

unid_summary <- unid %>% mutate(year_bin = floor(year / 5) * 5) %>% group_by(year_bin) %>% summarize(count_species = n())  # Count species records per bin

# Combine summaries and calculate proportions
combined_summary <- full_join(lithsum_summary, unid_summary, by = "year_bin") %>% replace_na(list(count_genus = 0, count_species = 0)) %>% mutate(total_count = count_genus,prop_genus = (count_genus - count_species) / total_count, prop_species = count_species / total_count)

#Convert to long format for ggplot
plot_data <- combined_summary %>% pivot_longer(cols = c(prop_genus, prop_species), names_to = "type",values_to = "proportion") %>% mutate(type = recode(type, prop_genus = "Specimens identified to genus", prop_species = "Specimens identified to species"))

# Plot with uniform bar height and labeled counts
gp <- ggplot(plot_data, ) +
  geom_bar(aes(x = year_bin, y = proportion, fill = type), stat = "identity", color = "white", width = 5) + #geom_area(aes(x = year_bin, y = total_count/7000),fill = "white",alpha = 0.8)+
  geom_line(aes(x = year_bin, y = total_count/9000), lty = 1,size=1.5, colour = "white")  + scale_fill_manual(values = c("Specimens identified to genus" = "#895C87", "Specimens identified to species" = "#B599B4")) +  theme_classic() + xlim(c(1850, 2020)) +
  xlab("Year of Collection") +  ylab("Proportion of Specimens") + theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.position = c(0.5, 1.05),
        legend.direction = "horizontal") +
  scale_y_continuous("Proportion of Specimens Indentified", sec.axis = sec_axis(~ . * 9000, name = "Number of Specimens Collected") 
  ) + theme(plot.margin =  unit(c(1.25,0.25,0,0.25),"cm"))

## save image file
## ggsave(gp , filename = "figure1_oct2024.png", units = "cm", height = 11,width = 17)
