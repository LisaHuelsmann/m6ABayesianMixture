

# Data prep: site level to transcript level


library(dplyr)




# Read data ---------------------------------------------------------------


dat_original = read.csv("data/UNCvsMETTL_DEGsTotal_m6A15Pct_PositionsAndRates.csv", sep = "\t")
head(dat_original)



# Calculate and explore site level methylation loss -----------------------

# calculate absolute and relative loss
dat_original %>% 
  filter(!is.na(MeanMeth_Unc)) %>%  # remove entries where MeanMeth_Unc is NA
  mutate(MeanMeth_loss = MeanMeth_Unc - MeanMeth_KD, # absolute loss
         MeanMeth_relloss = MeanMeth_loss/MeanMeth_Unc) ->  # relative loss
  dat

# viz
hist(dat$MeanMeth_loss)
hist(dat$MeanMeth_relloss, breaks = 40)
table(dat$MeanMeth_relloss >= 0.8)/nrow(dat)
table(dat$MeanMeth_relloss >= 1)/nrow(dat)
# most lose 100% anyway
# >= 80% relative loss would include 4% more sites

plot(log2FoldChange ~ MeanMeth_relloss, dat, ylim = c(-4, 4), 
     cex = 0.1, col = rgb(0, 0, 0, 0.01))
# no clear relationship between relative loss and log2FoldChange

# --> definition for "a site is lost": 100% methylation loss






# Aggregate transcript data -----------------------------------------------

dat %>% 
  group_by(gene_name) %>% 
  summarise(n_sites_all = n(), # number of sites
            n_sites_lost = sum(MeanMeth_relloss >= 1), # number of sites lost (define as 100% relative loss of methylation)
            mean_rel_location_all = mean(rel_location_total), # average location of sites
            mean_rel_location = mean(rel_location_total[MeanMeth_relloss >= 1]), # average location of lost sites
            # keep log2FC info
            log2FoldChange = unique(log2FoldChange),
            lfc_SE = unique(lfc_SE),
            log2FoldChange_shrunk = unique(log2FoldChange_shrunk),
            lfc_SE_shrunk = unique(lfc_SE_shrunk)) ->
  dat_transcript




# Save data ---------------------------------------------------------------

write.csv(dat_transcript, file = "data/Transcript_m6A_15.csv", row.names = F)
rm(list = ls())



