##### Jacana Microbiome Alpha Diversity Boxplots #####
library(tidyr)
library(dplyr)
require (officer)
require(rvg)
require(ggpubr)
library(ggplot2)

# Upload alpha diversity file

alphaData <- read.csv("alpha_metadata_jacana_micro.csv") #read in alpha diversity metadata
head(alphaData)
str(alphaData) #37 samples

# Change all independent variables to factors that appear in your metadata
alphaData$species <- as.factor(alphaData$species)
alphaData$sex <- as.factor(alphaData$sex)
alphaData$breeding_stage <- as.factor(alphaData$breeding_stage)

##### Plot for Model 1a: Species x Sex - FaithPD ##### 

p_faithpd_sex_species <- ggplot(alphaData, aes(x=species, y=FaithPD, fill=sex)) +
  geom_boxplot() + 
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab('Species') + 
  ylab('Faiths Phylogenetic Diversity (PD)') + 
  scale_x_discrete(labels=c('Wattled Jacana', 'Northern Jacana')) +
  scale_fill_manual(values=c('darksalmon', 'dodgerblue4'), 
                    labels = c('Female', 'Male'), name = "Sex") 

ggsave("faithpd_sex_species.pdf", height=4, width=5, device="pdf")


####### Plot dotplot Species x Sex - FaithPD
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

##### Plot for Model 1b: Species x Sex - Shannon ##### 

p_shannon_sex_species <- ggplot(alphaData, aes(x=species, y=Shannon, fill=sex)) +
  geom_boxplot() + 
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab('Species') + 
  ylab('Shannon Index') + 
  scale_x_discrete(labels=c('Wattled Jacana', 'Northern Jacana')) +
  scale_fill_manual(values=c('darksalmon', 'dodgerblue4'), 
                    labels = c('Female', 'Male'), name = "Sex") 

ggsave("shannon_sex_species.pdf", height=4, width=5, device="pdf")


####### Plot dotplot Species x Sex - Shannon
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

##### Plot for Model 1c: Species x Sex - Chao1 ##### 

p_chao1_sex_species <- ggplot(alphaData, aes(x=species, y=Chao1, fill=sex)) +
  geom_boxplot() + 
  theme_bw() +
  theme(text = element_text(size = 14), panel.grid.major = element_blank(),             panel.grid.minor = element_blank()) +
  xlab('Species') + 
  ylab('Chao1') + 
  scale_x_discrete(labels=c('Wattled Jacana', 'Northern Jacana')) +
  scale_fill_manual(values=c('darksalmon', 'dodgerblue4'), 
                    labels = c('Female', 'Male'), name = "Sex") 

ggsave("chao1_sex_species.pdf", height=4, width=5, device="pdf")

####### Plot dotplot Species x Sex - Chao1
p_chao1_sex_species_dot <- ggplot(alphaData, aes(x=species, y=Shannon, fill=sex)) +
  geom_boxplot() + geom_point(aes(shape=sex), position = position_jitterdodge(0.5)) +  
  ylab('Shannon Index') 
p_chao1_sex_species_dot

# Export plot to powerpoint
editable_graph <- dml(ggobj = p_chao1_sex_species_dot)
doc <- read_pptx()
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "p_chao1_sex_species_dot.pptx")

