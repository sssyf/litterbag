### microbial diversities
### Only codes for 16S are shown, the codes for ITS are similar.
library(ieggr)
library(picante)
library(lmerTest)


####################
### 1. read in data
####################
otu.file <- read.table("data/OTU_table_16S.txt", header = T, row.names = 1, check.names = F)
ITS.file <- read.table("data/OTU_table_ITS.txt", header = T, row.names = 1, check.names = F)
treat <- read.csv("data/treatment.csv")
rownames(treat) <- treat$sample

otu = t(otu.file)
ITS = t(ITS.file)



#################################
### 2. alpha diversity
#################################
### 2.1 16S
alpha_16S <- alpha.g(otu, trace = TRUE)
alpha_16S <- as.data.frame(alpha_16S)
alpha_16S$sample <- rownames(alpha_16S)
alpha_16S$shannon_diversity <- exp(1)^alpha_16S$shannon


### 2.2 linear mixed model
alpha_16S <- merge(treat, alpha_16S[, c(1,2,33)], by = "sample")
alpha_16S <- melt(alpha_16S, id = c("sample","year", "block", "Warm","type"))

alpha_16S_soil <- subset(alpha_16S, type == "soil")
alpha_16S_litter <- subset(alpha_16S, type == "litter")

alpha_name_16S <- unique(alpha_16S_soil$variable)
alpha_name_16S <- as.character(alpha_name_16S[c(1:9,34)])

# type
lm <- lmerTest::lmer(richness ~ type + (1|Warm) + (1|block) + (1|year), data = alpha_16S)
summary(lm)

# Warm - soil
lm <- lmerTest::lmer(richness ~ Warm + (1|block) + (1|year), data = subset(alpha_16S, type == "soil"))
summary(lm)

# Warm - litterbag
lm <- lmerTest::lmer(richness ~ Warm + (1|block) + (1|year), data = subset(alpha_16S, type == "litter"))
summary(lm)




#################################
### 3. beta diversity
#################################
### 3.1 PCoA
# Bray-Curtis
Bray_16S = vegdist(otu, method="bray")

pcoa.16S1 = pcoa(Bray_16S, correction = "none")
pcoa.16S1.samp = pcoa.16S1$vectors
pcoa.16S1.eig = pcoa.16S1$values
pcoa.16S1.sp = wascores(pcoa.16S1.samp[,1:3, drop=FALSE], otu, expand = TRUE)
pcoa.16S1.exp <- pcoa.16S1.eig$Relative_eig

# Sorensen
Soren_16S = vegdist(otu, method="bray", binary=TRUE)

pcoa.16S2 = pcoa(Soren_16S, correction = "none")
pcoa.16S2.samp = pcoa.16S2$vectors
pcoa.16S2.eig = pcoa.16S2$values
pcoa.16S2.sp = wascores(pcoa.16S2.samp[,1:3, drop=FALSE], otu, expand = TRUE)
pcoa.16S2.exp <- pcoa.16S2.eig$Relative_eig


### 3.2 dispersion
spc <- match.name(rn.list = list(otu = otu, treat = treat))
otu <- spc$otu
treat <- spc$treat
treat$group <- paste(treat$type, treat$Warm, sep = "_")

mod_16S1 <- betadisper(d = Bray_16S, group = treat$group, type = 'centroid')
mod_16S2 <- betadisper(d = Soren_16S, group = treat$group, type = 'centroid')

anova(mod_16S1)
TukeyHSD(mod_16S1, which = "group", ordered = FALSE, conf.level = 0.95)

anova(mod_16S2)
TukeyHSD(mod_16S2, which = "group", ordered = FALSE, conf.level = 0.95)



### 3.3 dissimilarity
# 3.3.1 soil/litter
i = 5
treat.use = treat[, i, drop=FALSE]
treat.name = colnames(treat)[i]

### 16S
dissim.16S.Soren = dissim(otu, treat = treat.use, dist.method = "bray", binary = TRUE)
colnames(dissim.16S.Soren) <- c("group1", "group2", "mrpp.delta", "mrpp.p", "anosim.R", "anosim.p", "adonis.F", "adonis.p")

dissim.16S.Bray = dissim(otu, treat = treat.use, dist.method = "bray", binary = FALSE)
colnames(dissim.16S.Bray) <- c("group1", "group2", "mrpp.delta", "mrpp.p", "anosim.R", "anosim.p", "adonis.F", "adonis.p")


# 3.3.2 soil - warming/unwarming
treat_soil <- subset(treat, type == "soil")
otu_soil <- subset(otu, rownames(otu) %in% treat_soil$sample)
otu_soil <- otu_soil[, colSums(otu_soil) != 0]

sampc = match.name(rn.list=list(otu_soil = otu_soil, treat_soil=treat_soil))
otu_soil = sampc$otu_soil
treat_soil = sampc$treat_soil

i = 4
treat.use_soil = treat_soil[, i, drop=FALSE]
treat.name_soil = colnames(treat_soil)[i]

dissim.16S.Soren_soil = dissim(otu_soil, treat = treat.use_soil, dist.method = "bray", binary = TRUE)
colnames(dissim.16S.Soren_soil) <- c("group1", "group2", "mrpp.delta", "mrpp.p", "anosim.R", "anosim.p", "adonis.F", "adonis.p")

dissim.16S.Bray_soil = dissim(otu_soil, treat = treat.use_soil, dist.method = "bray", binary = FALSE)
colnames(dissim.16S.Bray_soil) <- c("group1", "group2", "mrpp.delta", "mrpp.p", "anosim.R", "anosim.p", "adonis.F", "adonis.p")


# 3.3.3 litter - warming/unwarming
treat_litter <- subset(treat, type == "litter")
otu_litter <- subset(otu, rownames(otu) %in% treat_litter$sample)
otu_litter <- otu_litter[, colSums(otu_litter) != 0]

sampc = match.name(rn.list=list(otu_litter = otu_litter, treat_litter=treat_litter))
otu_litter = sampc$otu_litter
treat_litter = sampc$treat_litter

i = 4
treat.use_litter = treat_litter[, i, drop=FALSE]
treat.name_litter = colnames(treat_litter)[i]

dissim.16S.Soren_litter = dissim(otu_litter, treat = treat.use_litter, dist.method = "bray", binary = TRUE)
colnames(dissim.16S.Soren_litter) <- c("group1", "group2", "mrpp.delta", "mrpp.p", "anosim.R", "anosim.p", "adonis.F", "adonis.p")

dissim.16S.Bray_litter = dissim(otu_litter, treat = treat.use_litter, dist.method = "bray", binary = FALSE)
colnames(dissim.16S.Bray_litter) <- c("group1", "group2", "mrpp.delta", "mrpp.p", "anosim.R", "anosim.p", "adonis.F", "adonis.p")



### 3.4 Adonis
perm <- how(nperm = 999)
perm2 <- how(nperm = 999)
setBlocks(perm2) <- with(treat, block)

PERMANOVA_16S_Bray <- adonis2(otu ~ type*Warm*year, data = treat, method = "bray", permutations = perm2)
PERMANOVA_16S_Soren <- adonis2(otu ~ type*Warm*year, data = treat,  method = "bray", binary = TRUE, permutations = perm2)






