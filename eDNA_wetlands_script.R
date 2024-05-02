### Statistical analysis script in R v2023.06.1+524 

## Pearson’s product correlation tests 
# Data files are available on Wilderlab NZ explore map (https://www.wilderlab.co.nz/explore; accession codes 603456, 603161, 602697, 603227, 602353, 602509, 601982, 601772, 602406, 603499, 602000, 603082, 603498, 601570, 603118, 602725). Or, all data files are available in this Github repository. Note that the anonymous ('ANO') site was removed from all data files. 
#########################################
## Pearson’s product correlation tests ##
#########################################
## Total volume of water filtered
# Volume vs. Number of taxa
data <- read_excel("Correlation_Vol_Taxa_Count_WilderlabWetlandData.xlsx")
result <- cor.test(data$Volume, data$No_taxa)
print(result)
# Volume vs. Number of absolute read count
result2 <- cor.test(data$Volume, data$No_reads)
print(result2)
## Total number of replicates
# Replicates vs. Number of taxa
data2 <- read_excel("Correlation_Rep_Taxa_Count_WilderlabWetlandData.xlsx")
result3 <- cor.test(data2$Replicates, data2$No_taxa)
print(result3)
# Replicates vs. Number of absolute read count
result4 <- cor.test(data2$Replicates, data2$No_reads)
print(result4)

## Non-metric multidimensional plots
# Data files are available on Wilderlab NZ explore map (https://www.wilderlab.co.nz/explore; accession codes 603456, 603161, 602697, 603227, 602353, 602509, 601982, 601772, 602406, 603499, 602000, 603082, 603498, 601570, 603118, 602725). Or, all data files are available in this Github repository. Note that the anonymous ('ANO') site was removed from all data files.
#######################################
## Non-metric multidimensional plots ##
#######################################
## nMDS plot by site
# Load in the R packages
library(vegan)
library(pander)
library(ggplot2)
library(BiodiversityR)
library(readxl)
library(devtools)
# Read in the data
x1 <- read_excel("WetlandUIDs.xlsx", sheet = 3)
x1 <- x1[, c(1,8:ncol(x1))]
# Read in the metadata
metadata <- read_excel("WetlandUIDs.xlsx", sheet = 1)
# Working on merged data
community_matrix <- x1[seq(2, ncol(x1))]
UIDs <- colnames(community_matrix)
community_matrix <- t(as.matrix(community_matrix))
colnames(community_matrix) <- x1$Sequence
rownames(community_matrix) <- UIDs
colsums <- apply(community_matrix, 2, function(v) sum(v > 0))
community_matrix <- community_matrix[, colsums >= 7L]
mymds <- vegan::metaMDS(community_matrix, weakties = FALSE, autotransform = F)
mysamples <- as.data.frame(scores(mymds, "sites"))
mysamples$UID <- rownames(mysamples)
mysamples$Sites <- metadata$LocationCode[match(mysamples$UID, metadata$UID)]
myseqs <- as.data.frame(scores(mymds, "species"))
myseqs$species <- rownames(myseqs)
rownames(myseqs) <- NULL
# Order sites for ggplot 
mysamples$Sites <- factor(metadata$LocationCode, levels = c("ANO", "TAA", "HEN 1", "HEN 2", "MAT", "RAN", "OPU 1", "OPU 2", "OPU 3", "PON", "WKU", "NTT 1", "NTT 2", "NTT 3", "KAW 1", "KAW 2", "KAW 3", "NGA", "OLD", "QEP 1", "QEP 2", "INL", "OUT", "HOR", "RED", "RAK"))
# Plot sites with different colours
myplot <- ggplot(mysamples, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Sites), size = 2, stroke = 2.3) +
  scale_colour_manual(values = c("#922d4c", "#d83201", "#ad6a17", "#986432", "#6f9d54", "#6b8eb8", "#19a71f", "#a8e667", "#44f270", "#D6CBB5", "#CABEE9", "#926f61", "#a88141", "#d2b38c", "#2b4b3b", "#7D9D33", "#CED38C", "#f8cac2", "#daa7e7", "#c5d5f0", "#a1def0", "#f76015","#f79302", "#DCC949", "red", "#9f2108" )) + 
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(linewidth = 6),
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)
  ) +
  coord_equal()
myplot
## nMDS plot by location
# Replace site names with desired abbreviations in the data frame
mysamples$Sites <- gsub("HEN [12]", "HEN", mysamples$Sites)
mysamples$Sites <- gsub("OPU [123]", "OPU", mysamples$Sites)
mysamples$Sites <- gsub("NTT [123]", "NTT", mysamples$Sites)
mysamples$Sites <- gsub("KAW [123]", "KAW", mysamples$Sites)
mysamples$Sites <- gsub("QEP [12]", "QEP", mysamples$Sites)
# Order sites for ggplot
mysamples$Sites <- factor(mysamples$Sites, levels = c("ANO", "TAA", "HEN", "MAT", "RAN", "OPU", "PON", "WKU", "NTT", "KAW", "NGA", "OLD", "QEP", "INL", "OUT", "HOR", "RED", "RAK"))
# Plot locations with different colours
myplot <- ggplot(mysamples, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Sites), size = 2, stroke = 2.3) +
  scale_colour_manual(values = c("#922d4c", "#d83201", "#ad6a17", "#737a13", "#6b8eb8","#44f270", "#D6CBB5", "#CABEE9", "#926f61", "#2b4b3b", "#f8cac2", "#daa7e7", "#a1def0", "#f76015","#f79302", "#DCC949", "red", "#9f2108" )) + 
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(linewidth = 6),
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)
  ) +
  coord_equal()
myplot
## nMDS plot by region
# Assign sites into the three broad Aotearoa regions (i.e., Upper North Island, lower North Island, South Island)
metadata$Region <- ifelse(metadata$LocationCode %in% c("ANO", "TAA", "HEN 1", "HEN 2", "MAT", "RAN", "OPU 1", "OPU 2", "OPU 3", "PON", "WKU"), "Upper North Island",
                          ifelse(metadata$LocationCode %in% c("NTT 1", "NTT 2", "NTT 3", "KAW 1", "KAW 2", "KAW 3", "NGA", "OLD", "QEP 1", "QEP 2", "INL", "OUT"), "Lower North Island", 
                                 ifelse(metadata$LocationCode %in% c("HOR", "RED", "RAK"), "South Island", NA)))
# Order sites for ggplot
mysamples$Sites <- factor(metadata$Region, levels = c("Upper North Island", "Lower North Island", "South Island"))
# Plot regions with different colours
myplot <- ggplot(mysamples, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Sites), size = 2, stroke = 2.3) +
  scale_colour_manual(values = c("#6C90B9", "#313657", "#DD3C51")) +
  labs(color = "Region") +  # Set the legend title
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.key = element_rect(linewidth = 6),
    legend.key.size = unit(1, "cm"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14), 
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)
  ) +
  coord_equal()
myplot

## PERMANOVA and PERMDISP
# Data files are available on Wilderlab NZ explore map (https://www.wilderlab.co.nz/explore; accession codes 603456, 603161, 602697, 603227, 602353, 602509, 601982, 601772, 602406, 603499, 602000, 603082, 603498, 601570, 603118, 602725). Or, all data files are available in this Github repository. Note that the anonymous ('ANO') site was removed from all data files.
############################
## PERMANOVA and PERMDISP ##
############################
# Load in the R packages
library(vegan)
# Formatting frequency table
pres_abs <- read.csv("Wilderlab_perm.csv")
pres_abs1 <- pres_abs[,1] 
pres_abs2 <- read.csv("Wilderlab_perm.csv")
pres_abs3 <- pres_abs2[,2:27]
pres_abs3[pres_abs3 > 0] <- 1
pres_abs_final <- cbind(pres_abs1, pres_abs3)
perm <- t(pres_abs_final) 
colnames(perm) <- perm[1,]
perm <- perm[-1, ]
write.csv(perm, file="finalmat_perm.csv", row.names=T)
# Read in the data
metadata <- read.csv("Wilderlab_perm_meta_ULS.csv")
permanova <- read.csv("finalmat_perm.csv", header = TRUE, row.names = 1)
# Run PERMANOVA
adonis2(permanova ~ Cluster, data = metadata, by = 'margin')
# Test for dispersion effects
dis <- vegdist(permanova, method = 'jaccard')
groups <- factor(c(rep(1,10), rep(2,13), rep(3,3)),
                 labels = c('Upper', 'Lower', 'South'))
dis
groups
mod <- betadisper(dis, groups)
mod
boxplot(mod)
set.seed(25)
permutest(mod)
plot(mod)
set.seed(4)
permutest(mod, pairwise = TRUE)
