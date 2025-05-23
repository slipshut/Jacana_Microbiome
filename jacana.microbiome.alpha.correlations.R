# Jacana phenotypic integration with microbial diversity
library(tidyr)
library(dplyr)
require (officer)
require(rvg)
require(ggpubr)
library(ggplot2)

# Set working directory
setwd("~/Dropbox/Duke/Research Projects/Jacana microbiome")

# Read in csv file with columns for species, sex, individual, and traits
jacana <- read.csv("metadata_jacana_micro.csv", header=TRUE, sep=",")

# Visualize data
View(jacana)
ggplot(aes(y=Chao1, x=species), data = jacana) + geom_boxplot() + scale_fill_brewer(palette = "Reds") + theme_classic()

# Subset groups by species and sex
wattled.female = subset (jacana, species == "J. jacana" & sex == "female")
wattled.male = subset (jacana, species == "J. jacana" & sex == "male")
northern.female = subset (jacana, species == "J. spinosa" & sex == "female")
northern.male = subset (jacana, species == "J. spinosa" & sex == "male")

# Subset groups by traits 
wattled.female.int <- wattled.female[, c("log_testosterone","avg_spur","body_mass", "Chao1", "Shannon","FaithPD")]
wattled.male.int <- wattled.male[, c("log_testosterone","avg_spur","body_mass", "Chao1", "Shannon","FaithPD")]
northern.female.int <- northern.female[, c("log_testosterone","avg_spur","body_mass", "Chao1", "Shannon","FaithPD")]
northern.male.int <- northern.male[, c("log_testosterone","avg_spur","body_mass", "Chao1", "Shannon","FaithPD")]

# Use complete observations for correlations
res.wattled.female <- cor(wattled.female.int, method = "spearman", use = "complete.obs")
res.wattled.male <- cor(wattled.male.int, method = "spearman", use = "complete.obs")
res.northern.female <- cor(northern.female.int, method = "spearman", use = "complete.obs")
res.northern.male <- cor(northern.male.int, method = "spearman", use = "complete.obs")

# # visualize correlations, including nonsignificant ones
# library(corrplot)
# 
# corrplot(res.wattled.female, type = "upper", order = "hclust", )
# corrplot(res.wattled.male, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
# corrplot(res.northern.female, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
# corrplot(res.northern.male, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)


# # mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(wattled.female.int)
corrplot(res.wattled.female, type="upper", p.mat = p.mat, tl.col = "black", tl.srt = 45, sig.level = 0.05, insig = "blank")

p.mat <- cor.mtest(wattled.male.int)
corrplot(res.wattled.male, type="upper", p.mat = p.mat, tl.col = "black", tl.srt = 45, sig.level = 0.05, insig = "blank")

p.mat <- cor.mtest(northern.female.int)
corrplot(res.northern.female, type="upper", p.mat = p.mat, tl.col = "black", tl.srt = 45, sig.level = 0.05, insig = "blank")

p.mat <- cor.mtest(northern.male.int)
corrplot(res.northern.male, type="upper", p.mat = p.mat, tl.col = "black", tl.srt = 45, sig.level = 0.05, insig = "blank")

# Visualize raw data trends
# northern jacana females
northern.female = subset (jacana, species == "J. spinosa" & sex == "female")

ggplot(northern.female, aes(x=Chao1, y=distance)) + geom_point( color="#69b3a2")

ggplot(northern.female, aes(x=log_testosterone, y=avg_spur)) + geom_point( color="#69b3a2")
cor.test(northern.female$log_testosterone, northern.female$avg_spur) #t = 3.3554, df = 8, p-value = 0.01, cor = 0.7645929 

ggplot(northern.female, aes(x=Chao1, y=avg_spur)) + geom_point( color="#69b3a2")
cor.test(northern.female$Chao1, northern.female$avg_spur) #t = 1.3487, df = 9, p-value = 0.2104, cor = 0.4100312 

ggplot(northern.female, aes(x=log_testosterone, y=Chao1)) + geom_point(size = 2.5) + 
  theme_classic() + geom_smooth(method=lm, color = "black") 
cor.test(northern.female$log_testosterone, northern.female$Chao1, method= "pearson") #t = 3.131, df = 8, p-value = 0.01399, cor = 0.7420518 

ggplot(northern.female, aes(x=log_testosterone, y=FaithPD)) + geom_point(size = 2.5) + 
  theme_classic() + geom_smooth(method=lm, color = "black") 
cor.test(northern.female$log_testosterone, northern.female$FaithPD) #t = 2.5087, df = 8, p-value = 0.03644, cor = 0.6635639 

ggplot(northern.female, aes(x=FaithPD, y=avg_spur)) + geom_point( color="#69b3a2")
cor.test(northern.female$FaithPD, northern.female$avg_spur) #t = 1.1814, df = 9, p-value = 0.2677, cor = 0.3664164


# northern jacana males
northern.male = subset (jacana, species == "J. spinosa" & sex == "male")
ggplot(northern.male, aes(x=log_testosterone, y=body_mass)) + geom_point( color="#69b3a2")
cor.test(northern.male$log_testosterone, northern.male$body_mass) #t = 3.022, df = 10, p-value = 0.01285, cor = 0.6908852

ggplot(northern.male, aes(x=body_mass, y=Shannon, color=breeding_stage, shape=breeding_stage)) + geom_point(size = 2.5) + 
  geom_smooth(method=lm, aes(fill=breeding_stage)) + theme_classic() 
cor.test(northern.male$body_mass, northern.male$Shannon) #t = 2.4205, df = 10, p-value = 0.03603, cor = 0.6078077 

ggplot(northern.male, aes(x=log_testosterone, y=Shannon, color=breeding_stage, shape=breeding_stage)) + geom_point()
cor.test(northern.male$log_testosterone, northern.male$Shannon) #t = 0.48058, df = 10, p-value = 0.6412, cor = 0.1502483 

