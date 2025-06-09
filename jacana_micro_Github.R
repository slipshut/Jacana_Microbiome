### Jacana Microbiome Analysis ###
# Last updated: 6/9/2025

# Load libraries

# pacman allows you to install multiple packages at the same time rather than have to use library() for each package individually

pacman::p_load("phyloseq", "ggplot2", "viridis", "grid", "gridExtra", "lme4", "lmerTest", "emmeans","vegan", "devtools", "sjPlot", "here", "ape",
               "ggordiplots", "BiodiversityR", "dplyr", "tidyr", "DHARMa", "microbiome")

#Upload ASV, taxonomy, tree and metadata files to phyloseq:

#ASV table
unrarefied_ASV_table <- read.csv("unrarefied_ASV_table.csv", row.names=1, sep = ",") 

# read in asv table with feature names as rownames
# str() gives summary of the data
str (unrarefied_ASV_table) #2441 obs. of 43 variables, revised 2231 obs. of 43
unrarefied_ASV_table <- as.matrix (unrarefied_ASV_table) 
#make into a matrix
View(unrarefied_ASV_table)

#taxonomy file

taxonomy <- read.csv ("taxonomy.csv", row.names = 1)
taxonomy <- separate(taxonomy,Taxon,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep = "; ")
str(taxonomy) # 2445 obs.
taxonomy <- as.matrix (taxonomy)

#read in tree as a phyloseq object

phy_tree <- read.tree("tree.nwk")

#load metadata file

metadatafull <- read.csv("metadata_jacana_micro.csv")
str(metadatafull)  # contains 43 samples (has to match ASV table observations)
row.names(metadatafull)<-metadatafull$sample_name

#import all as phyloseq objects
ASV <- otu_table(unrarefied_ASV_table, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadatafull)
head(META) # head() gives the first 5 lines of the file to glance over

#check that the ASV and sample names are consistent across objects

str(taxa_names(TAX))
str(taxa_names(ASV))
str(taxa_names(phy_tree))

str(sample_names(ASV))
str(sample_names(META))

#### MERGE INTO PHYLOSEQ OBJECT ####
physeq <- phyloseq(ASV, TAX, META, phy_tree) # physeq() is taking the ASV, TAX, META, phy_tree and combining to a physeq object
physeq #2231 taxa, 43 samples 

#check the number of reads per sample
sample_reads<-data.frame(reads=sample_sums(physeq))
View(sample_reads)
write.csv(sample_reads, "sample_read_counts.csv")
# Samples with lowest reads: WAF3 (7445), WAM6 (221)

##### Decontam #####

# Identifying which bacterial DNA is from negative kit controls (Blanks) vs an actual sample (true_sample)

library(decontam)

# Remove contaminants with the combined function of prevalence and frequency
sample_data(physeq)$is.neg <- sample_data(physeq)$sample_or_control == "negative_control"
contamdf.prev <- isContaminant(physeq, method="combined",conc="quant_reading",neg="is.neg")
table(contamdf.prev$contaminant)

# Creating the decontaminated phyloseq object (physeq_nc)
physeq_nc<-prune_taxa(!contamdf.prev$contaminant,physeq) #physeq_nc is how much after removing ASVs found in blanks

physeq_nc #identified 4 contaminants and removed them leaving 2227 taxa

# Remove negative controls (n=4)

# Removing the Blanks
physeq2 <- subset_samples(physeq_nc, sample_or_control == "true_sample") 
physeq2 #39 samples, 2227 taxa

##### RAREFY READS TO MIN SAMPLING DEPTH ######

# Look at lowest sample depth
sample_sums(physeq2)

#rarefy to 14000 and set seed before rarefying (28367), so that results are reproducible  

physeqRare<-rarefy_even_depth(physeq2, 14000, rngseed = 28367)
# We do not want samples with such low reads because it can throw our data off

# leaves 1989 taxa across 37 samples
physeqRare
sample_sums(physeqRare)

##### CALCULATE ALPHA DIVERSITY METRICS ####

# calculate chao1, shannon using estimate_richness()
# richness estimates for each sample, with measures of Chao1 and Shannon
set.seed(787)
richnessEstRare<-estimate_richness(physeqRare, split=TRUE, measures= c("Chao1", "Shannon"))
head(richnessEstRare)
str(richnessEstRare)

# Calculating Faith's PD with picante package

library(picante)

faiths_asv <- as.data.frame(physeqRare@otu_table) # creates dataframe from otu table
faiths_tree <- physeqRare@phy_tree # extracts phylogenetic tree into faiths_tree
faithPD <- pd(t(faiths_asv),faiths_tree, include.root = T) # t transposes the otu table for use in picante; include.root = T(rue) includes root node in calculation
head(faithPD)
richnessEstRare$FaithPD <- faithPD$PD 
head(richnessEstRare)
write.csv(richnessEstRare, "richnessEstRare.csv")

#add alpha diversity metrics to metadata

physeqRareMeta <- as.data.frame(sample_data(physeqRare))
head(physeqRareMeta)
physeqRareMeta$Shannon <- richnessEstRare$Shannon
physeqRareMeta$Chao1 <- richnessEstRare$Chao1
physeqRareMeta$FaithPD <- richnessEstRare$FaithPD
head(physeqRareMeta)
str(physeqRareMeta)

alphaData <- read.csv("metadata_jacana_micro.csv") #read in metadata that has alpha diversity metrics added
head(alphaData)
str(alphaData)

# Exclude the blanks and samples removed during rarefaction (WAF3, WAM6)
alphaData <-subset(alphaData,sample_or_control!="negative_control") 
alphaData <-subset(alphaData,sample_name!="WAF3") 
alphaData <-subset(alphaData,sample_name!="WAM6") 

# Change all independent variables to factors that appear in your metadata
alphaData$species <- as.factor(alphaData$species)
alphaData$sex <- as.factor(alphaData$sex)

##### Alpha Diversity Figures #####
library("officer")
library("rvg")
library("ggpubr")

##### Figure S1A Species x Sex - Chao1 ##### 

p_chao1_sex_species_dot <- ggplot(alphaData, aes(x=species, y=Chao1, fill=sex)) +
  geom_boxplot() + geom_point(aes(shape=sex), position = position_jitterdodge(0.5)) +  
  ylab('Chao1 Index') 
p_chao1_sex_species_dot

# Export plot to powerpoint
editable_graph <- dml(ggobj = p_chao1_sex_species_dot)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "p_chao1_sex_species_dot.pptx")

##### Figure S1B Species x Sex - Shannon ##### 

p_shannon_sex_species_dot <- ggplot(alphaData, aes(x=species, y=Shannon, fill=sex)) +
  geom_boxplot() + geom_point(aes(shape=sex), position = position_jitterdodge(0.5)) +  
  ylab('Shannon Index') 
p_shannon_sex_species_dot

# Export plot to powerpoint
editable_graph <- dml(ggobj = p_shannon_sex_species_dot)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "p_shannon_sex_species_dot.pptx")

##### Figure S1C Species x Sex - FaithPD ##### 

p_faithpd_sex_species_dot <- ggplot(alphaData, aes(x=species, y=FaithPD, fill=sex)) +
  geom_boxplot() + geom_point(aes(shape=sex), position = position_jitterdodge(0.5)) +  
  ylab('Faiths Phylogenetic Diversity (PD)') 
p_faithpd_sex_species_dot

# Export plot to powerpoint
editable_graph <- dml(ggobj = p_faithpd_sex_species_dot)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "p_faithpd_sex_species_dot.pptx")

##### Alpha Diversity Interaction Models #####

##### Table S1 LMM Species x Sex - Chao1 #####

mod_chao1_sex_species <- lm(Chao1 ~ species*sex, data = alphaData)
summary(mod_chao1_sex_species, ddf = "Kenward-Roger")
tab_model(mod_chao1_sex_species, p.val = "kr", show.df = TRUE, file = "mod_results_chao1_sex_species.html")

#model diagnostics with DHARMa

simulationOutput_chao1_sex_species <- simulateResiduals(fittedModel = mod_chao1_sex_species, plot = F)
plot(simulationOutput_chao1_sex_species)
plotResiduals(simulationOutput_chao1_sex_species, form = alphaData$Chao1)
testDispersion(simulationOutput_chao1_sex_species)

##### Table S2 LMM Species x Sex - Shannon #####

mod_shannon_sex_species <- lm(Shannon ~ species*sex, data = alphaData)
summary(mod_shannon_sex_species, ddf = "Kenward-Roger")
tab_model(mod_shannon_sex_species, p.val = "kr", show.df = TRUE, file = "mod_results_shannon_sex_species.html")

#model diagnostics with DHARMa

simulationOutput_shannon_sex_species <- simulateResiduals(fittedModel = mod_shannon_sex_species, plot = F)
plot(simulationOutput_shannon_sex_species)
plotResiduals(simulationOutput_shannon_sex_species, form = alphaData$Shannon)
testDispersion(simulationOutput_shannon_sex_species)

##### Table S3 LMM Species x Sex - FaithPD #####

mod_faithpd_sex_species <- lm(FaithPD ~ species*sex, data = alphaData)
summary(mod_faithpd_sex_species, ddf = "Kenward-Roger")
tab_model(mod_faithpd_sex_species, p.val = "kr", show.df = TRUE, file = "mod_results_faithpd_sex_species.html")

#model diagnostics with DHARMa

simulationOutput_faithpd_sex_species <- simulateResiduals(fittedModel = mod_faithpd_sex_species, plot = F)
plot(simulationOutput_faithpd_sex_species)
plotResiduals(simulationOutput_faithpd_sex_species, form = alphaData$FaithPD)
testDispersion(simulationOutput_faithpd_sex_species)

##### BETA DIVERSITY #####

##### PERMANOVA - Bray Curtis sex*species #####

# Use the rarefied phyloseq object physeqRare

### RAREFIED PRUNE RARE TAXA ######

# Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
prevdf_rare = apply(X = otu_table(physeqRare),
                    MARGIN = ifelse(taxa_are_rows(physeqRare), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts for each phylum to this data.frame
prevdf_rare = data.frame(Prevalence = prevdf_rare,
                         TotalAbundance = taxa_sums(physeqRare),
                         tax_table(physeqRare))
head(prevdf_rare)
str(prevdf_rare)

plyr::ddply(prevdf_rare, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum

# Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
prevdf1_rare = subset(prevdf_rare, Phylum %in% get_taxa_unique(physeqRare, "Phylum"))
ggplot(prevdf1_rare, aes(TotalAbundance, Prevalence / nsamples(physeqRare),color=Phylum)) +
  # Include filtering threshold here; 2%, total abundance across all samples set to 50
  geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define abundance as 50 total reads across samples

abundanceThreshold<-50

# Execute the prevalence filter, using `prune_taxa()` function

head(prevdf1_rare)

KeepTaxa1_rare<- rownames(prevdf1_rare)[(prevdf1_rare$TotalAbundance >= abundanceThreshold)]
str(KeepTaxa1_rare) #352 taxa

physeqBeta_rare<- prune_taxa(KeepTaxa1_rare,physeqRare)
physeqBeta_rare #352 taxa

##### Export rarefied ASV table and taxonomy into qiime2 for taxabarplot Figure S2 #####

#Export taxonomy table as "tax.txt"

tax<-as(tax_table(physeqBeta_rare),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "tax_taxa_barplot.txt", quote=FALSE, col.names=FALSE, sep="\t")

#Export feature/OTU table as a biom table
######DO NOT TRANSFORM OR THE FEATURES AND SAMPLES WILL BE SWITCHED

library(biomformat);packageVersion("biomformat")

otu<-as(otu_table(physeqBeta_rare),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_taxa_barplot.biom")

#Extract metadata

beta_data<- data.frame(sample_data(physeqBeta_rare))

# Look at structure and make sure they are correct data types

str(beta_data)
beta_data$species <- as.factor(beta_data$species)
beta_data$sex <- as.factor(beta_data$sex)

##### Bray-Curtis Beta Diversity #####
##### Figure S2 PCoA Bray - sex and species#####

ordu_bray = ordinate(physeqBeta_rare, "PCoA", "bray", weighted=FALSE)
P_bray <- plot_ordination(physeqBeta_rare, ordu_bray, color="species", shape="sex")
pca_sex_species_bray <- P_bray + geom_point(alpha=0.5, size=3) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(values=c(16, 17), labels = c("Female", "Male"), name = "Sex") +
  scale_color_manual(values=c('red2', 'goldenrod1'), labels = c("J.jacana", "J.spinosa"), name = "Species")
pca_sex_species_bray
ggsave("pca_sex_species_bray.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

plot_ordination(physeqBeta_rare, ordu_bray, color="species", shape="sex")

# Export plot to powerpoint
editable_graph <- dml(ggobj = pca_sex_species_bray)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "pca_sex_species_bray.pptx")

# Calculate Bray-Curtis dissimiliary
bray <- phyloseq::distance(physeqBeta_rare, method = "bray")

bray_matrix <- as.matrix(dist(bray))

write.csv(bray_matrix, "bray_matrix.csv")

##### PERMANOVA Bray - Sex and Species #####

perm <- how(nperm = 9999)
set.seed(5498)
permanova_sex_species_bray<- adonis2(bray ~ sex*species, data=beta_data, 
                                   permutations = perm, method = "bray", by = "margin")
permanova_sex_species_bray

#Export results

str(permanova_sex_species_bray)
perm_results_sex_species_bray <- data.frame(permanova_sex_species_bray[c(1,2,3,4,5)])
perm_results_sex_species_bray
write.csv(perm_results_sex_species_bray, "permanova_results_sex_species_bray.csv")

# Sex only 

perm <- how(nperm = 9999)
set.seed(5498)
permanova_sex_bray<- adonis2(bray ~ sex, data=beta_data, 
                                     permutations = perm, method = "bray", by = "margin")
permanova_sex_bray

#Export results

str(permanova_sex_bray)
perm_results_sex_bray <- data.frame(permanova_sex_bray[c(1,2,3,4,5)])
perm_results_sex_bray
write.csv(perm_results_sex_bray, "permanova_results_sex_bray.csv")

# Species only 

perm <- how(nperm = 9999)
set.seed(5498)
permanova_species_bray<- adonis2(bray ~ species, data=beta_data, 
                             permutations = perm, method = "bray", by = "margin")
permanova_species_bray

#Export results

str(permanova_species_bray)
perm_results_species_bray <- data.frame(permanova_species_bray[c(1,2,3,4,5)])
perm_results_species_bray
write.csv(perm_results_species_bray, "permanova_results_species_bray.csv")

##### PERMDISP Sex and Species #####
### Checking group dispersions with betadisper analysis

# Sex

beta_disp_sex_bray<-betadisper(bray, beta_data$sex)
set.seed(25)
permdisp_sex_bray<- permutest(beta_disp_sex_bray, permutations = 9999)
permdisp_sex_bray

# Species

beta_disp_species_bray<-betadisper(bray, beta_data$species)
set.seed(25)
permdisp_species_bray<- permutest(beta_disp_species_bray, permutations = 9999)
permdisp_species_bray

##### Jaccard Beta Diversity #####
##### Figure 2 PCoA Jaccard - sex and species#####

ordu_jaccard = ordinate(physeqBeta_rare, "PCoA", "jaccard", weighted=FALSE)
P_jaccard <- plot_ordination(physeqBeta_rare, ordu_jaccard, color="species", shape="sex")
pca_sex_species_jaccard <- P_jaccard + geom_point(alpha=0.5, size=3) + #alpha controls transparency and helps when points are overlapping
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_shape_manual(values=c(16, 17), labels = c("Female", "Male"), name = "Sex") +
  scale_color_manual(values=c('red2', 'goldenrod1'), labels = c("J. jacana", "J. spinosa"), name = "Species")
pca_sex_species_jaccard
ggsave("pca_sex_species_jaccard.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

plot_ordination(physeqBeta_rare, ordu_jaccard, color="species", shape="sex")


# Export plot to powerpoint
editable_graph <- dml(ggobj = pca_sex_species_jaccard)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "pca_sex_species_jaccard.pptx")

# Calculate Jaccard dissimiliary
jaccard <- phyloseq::distance(physeqBeta_rare, method = "jaccard")

jaccard_matrix <- as.matrix(dist(jaccard))

write.csv(jaccard_matrix, "jaccard_matrix.csv")

##### PERMANOVA Jaccard - Sex and Species #####

perm <- how(nperm = 9999)
set.seed(5498)
permanova_sex_species_jaccard<- adonis2(jaccard ~ sex*species, data=beta_data, 
                                     permutations = perm, method = "jaccard", by = "margin")
permanova_sex_species_jaccard

#Export results

str(permanova_sex_species_jaccard)
perm_results_sex_species_jaccard <- data.frame(permanova_sex_species_jaccard[c(1,2,3,4,5)])
perm_results_sex_species_jaccard
write.csv(perm_results_sex_species_jaccard, "permanova_results_sex_species_jaccard.csv")

# Sex only 

perm <- how(nperm = 9999)
set.seed(5498)
permanova_sex_jaccard<- adonis2(jaccard ~ sex, data=beta_data, 
                             permutations = perm, method = "jaccard", by = "margin")
permanova_sex_jaccard

#Export results

str(permanova_sex_jaccard)
perm_results_sex_jaccard <- data.frame(permanova_sex_jaccard[c(1,2,3,4,5)])
perm_results_sex_jaccard
write.csv(perm_results_sex_jaccard, "permanova_results_sex_jaccard.csv")

# Species only 

perm <- how(nperm = 9999)
set.seed(5498)
permanova_species_jaccard<- adonis2(jaccard ~ species, data=beta_data, 
                                 permutations = perm, method = "jaccard", by = "margin")
permanova_species_jaccard

#Export results

str(permanova_species_jaccard)
perm_results_species_jaccard <- data.frame(permanova_species_jaccard[c(1,2,3,4,5)])
perm_results_species_jaccard
write.csv(perm_results_species_jaccard, "permanova_results_species_jaccard.csv")

##### PERMDISP Sex and Species - Jaccard #####
### Checking group dispersions with betadisper analysis

# Sex

beta_disp_sex_jaccard<-betadisper(jaccard, beta_data$sex)
set.seed(25)
permdisp_sex_jaccard<- permutest(beta_disp_sex_jaccard, permutations = 9999)
permdisp_sex_jaccard

# Species

beta_disp_species_jaccard<-betadisper(jaccard, beta_data$species)
set.seed(25)
permdisp_species_jaccard<- permutest(beta_disp_species_jaccard, permutations = 9999)
permdisp_species_jaccard

##### Table S6 Maaslin2 Analysis - Bacteria + Testosterone in Female Northern Jacanas #####
# Need to match genus codes from the taxonomy file

library("Maaslin2")

#filter beta diversity to only include females
female_physeqBeta_rare<- subset_samples(physeqBeta_rare,sex =="female")

#filter beta diversity to only include northern jacanas
NO_female_physeqBeta_rare<- subset_samples(female_physeqBeta_rare,species=="J. spinosa")

#filter beta diversity to remove NOF6 which doesn't have a T measurement
NO_female_physeqBeta_rare<- subset_samples(NO_female_physeqBeta_rare, !is.na(log_testosterone))

# Agglomorate to genus level
NO_female_physeqBeta_rare <- tax_glom(NO_female_physeqBeta_rare, "Genus")
NO_female_physeqBeta_rare

# Extract the taxonomy table from the filtered phyloseq object to find taxa names
Maaslin2_tax_table <- tax_table(NO_female_physeqBeta_rare)
write.csv(Maaslin2_tax_table, "Maaslin2_tax_table.csv")

# Run Maaslin2 to correlate bacterial genus relative abundance with log testosterone levels in female northern jacanas
mas_1 <- Maaslin2(
  input_data = data.frame(otu_table(NO_female_physeqBeta_rare, taxa_are_rows = TRUE)),
  input_metadata = data.frame(sample_data(NO_female_physeqBeta_rare)),
  output = "./Maaslin2_female_northern_jacana_analysis_experimental",
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "TSS",
  transform = "AST",
  analysis_method = "CPLM",
  max_significance = 0.05,
  fixed_effects = "log_testosterone",
  correction = "BH",
  standardize = FALSE,
  cores = 1,
  plot_heatmap = TRUE,
  heatmap_first_n = 20,
  plot_scatter = TRUE,
  max_pngs = 10,
  save_scatter = FALSE,
  save_models = FALSE,
  reference = NULL)
mas_res_df <- mas_1$results

fdr_mas <- mas_res_df %>%
  dplyr::filter(qval < 0.05)

dim(fdr_mas)
head(fdr_mas)
