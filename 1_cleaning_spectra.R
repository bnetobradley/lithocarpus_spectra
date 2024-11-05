#### READING AND CLEANING SPECTRA ####

## READ IN LIBRARIES
library(spectrolab)
library(tidyverse)

## read in metadata
met.df <- read.csv("data/04032022_lithocarpus/working_metadata-Table 1.csv")
met.df$samplename <- paste0(met.df$spectra_file, ".sig")

## read in spectral data
specdat <- read_spectra(path = "data/spectra_from_kew/final_lithocarpus_spectra/", format = "sig")

## remove reflectance values greater than 1
specdat <- specdat[!rowSums(specdat > 1),]

## match the data around the sensors
specdat <- match_sensors(specdat, splice_at = c(1000,1830))

## cut out wavelengths beyond 400 and 2400
wave <- wavelengths(specdat)
wave <- wave[wave < 2400 & wave > 400]
specdat <- specdat[, wave]

## smooth the spectra
spectdat <- smooth(specdat)

## join spectra with metadata
spec.df <- as.data.frame(specdat)
spec.df <- inner_join(spec.df, met.df, by = c("sample_name" = "samplename"))

## cut out scans noted as bad
spec.df <- spec.df %>% filter(leaf_description != "BAD SCAN" & leaf_description != "bad scan" & leaf_description != "bAD SCAN")

## remove non-lithocarpus scans
spec.df <- spec.df %>% filter(Genus != "Maesa" & Genus != "Erycibe") %>% filter(species != "glomerata" & species != " macrocarpa" & species != " bismarckiana" & species != "grandifolia")


## turn back to spectra
metadat <- as.data.frame(spec.df[,904:922])
metadat$samplename <- spec.df$sample_name
specdat <- as.spectra(spec.df[,1:903], name_idx = 1)

## save as .rds
#saveRDS(specdat, "data/clean_lithocarpus_spectra.rds")
#saveRDS(metadat, "data/clean_lithocapurs_metadat.rds")

## vector normalize spectra
specdat <- normalize(specdat)

#saveRDS(specdat, "data/clean_normalized_lithocarpus_spectra.rds")