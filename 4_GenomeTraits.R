### genome traits
### Only codes for 16S are shown, the codes for ITS are similar.
library(readr)
library(dplyr)
library(ieggr)


### 1. download prokaryotic genome assemblies from NCBI, and do blast in our server
### 2. extract the best hit (besthit.txt)

#############################################
### 3. calculate genome traits for each zOTU
#############################################
besthit <- read_delim("data/besthit.txt")
besthit <- besthit[, 1:7]

# some sequences were partially hit to whole genome, with low bit.score but high identity, here is to filter those hits 
# extract the highest hit.score of each sequence
besthit1 <- subset(besthit, bit.score >= 100 & q.length >= 100)
length(unique(besthit1$query.id))

summary(besthit1$bit.score)
summary(besthit1$q.length)
summary(besthit1$e.value)
summary(besthit1$identity)


# seq id with max hit.score and 100%, 99.5%, 98.7%, and 94.5% identity
m <- subset(besthit1, select=c('query.id','bit.score'))
# seq id and max hit.score
m1 <- aggregate(. ~ query.id, m, max)
# all information of the seq id with max hit.score
besthit2 <- merge(besthit1, m1, by=c('query.id','bit.score'), all.y=T, sort=F)
m100 <- unique(subset(besthit2, identity==100 & bit.score>=400, select=c(query.id)))   #781
m99.5 <- unique(subset(besthit2, identity>=99.5 & bit.score>=400, select=c(query.id))) #1245
m98.7 <- unique(subset(besthit2, identity>=98.7 & bit.score>=400, select=c(query.id))) #2478
m94.5 <- unique(subset(besthit2, identity>=94.5 & bit.score>=400, select=c(query.id))) #12124


# rrn copy number
besthit3 <- besthit1
besthit3$bit.score <- besthit3$e.value <- besthit3$identity <- besthit3$q.length <- NULL
besthit3$hit <- rep(1,nrow(besthit3))

besthit4 <- subset(besthit3, select=c('query.id','Assembly.Accession','hit'))
besthit4 <- aggregate(. ~ query.id + Assembly.Accession, besthit4, sum)
str(besthit4)

besthit5 <- aggregate(hit ~ query.id, besthit4, mean)
summary(besthit5$hit)


# genome size
prokaryotes <- read.table("data/prokaryotes2.txt", header = T, sep = ";")
colnames(prokaryotes)
gsize <- merge(besthit4, prokaryotes[, c(1,5:8,12:13,19)], by = "Assembly.Accession")
str(gsize)

colnames(gsize)
genomesize <- gsize[, c(1,2,7,8)]
genomesize2 <- aggregate(Size..Mb. ~ query.id, genomesize, mean)

gsize$GC. <- as.numeric(as.character(gsize$GC.))
gsize2 <- subset(gsize, is.na(gsize$GC.)==F)
for (i in genomesize2$query.id) {
  gsizei <- subset(gsize2, query.id == i)
  genomesize2[genomesize2$query.id == i, "GC"] <- sum(gsizei$Size..Mb. * gsizei$GC.) / sum(gsizei$Size..Mb.)
} 
summary(genomesize2$GC)
colnames(genomesize2)[2] <- "size"

genometrait <- merge(besthit5, genomesize2, by = "query.id")
colnames(genometrait)[1] <- "OTU.ID"




################################################
### 4. calculate community-level genome traits
################################################
### 4.1 read in data
treat <- read.csv("data/treatment.csv")
rownames(treat) <- treat$sample
treat$group <- paste(treat$type, treat$Warm, sep = "_")

otu.file <- read.table("data/OTU_table_16S.txt", header = T, row.names = 1, check.names = F)
otu = as.data.frame(t(otu.file))


### 4.2 calculate
samp <- match.name(rn.list = list(otu = otu, treat = treat))
otu <- samp$otu
treat <- samp$treat

otu <- as.data.frame(t(otu))
otu$OTU.ID <- rownames(otu)
otu <- merge(otu, genometrait, by = "OTU.ID")

trait <- as.data.frame(colnames(otu)[2:(ncol(otu)-3)])
colnames(trait) <- "sample"
for (i in 1:48) {
  trait[i, "CopyNumber"] <- sum(otu[, i+1]) / sum(otu[, i+1]/otu[, 50])
}

for (i in 1:48) {
  trait[i, "GenomeSize"] <- sum(otu[, i+1]*otu[, 51]) / sum(otu[, i+1])
}

for (i in 1:48) {
  trait[i, "GC"] <- sum(otu[, i+1]*otu[, 52]) / sum(otu[, i+1])
}


### 4.3 compare
trait <- merge(treat, trait, by = "sample")

## copy number
anova(lm(CopyNumber~group, trait))
trait_aov <- aov(lm(CopyNumber~group, trait))
TukeyHSD(trait_aov)
anova(lm(CopyNumber~type, trait))
t.test(subset(trait, type == "litter")$CopyNumber, subset(trait, type == "soil")$CopyNumber)

summary(lmerTest::lmer(CopyNumber ~ type * Warm + (1|block) + (1|year), data = trait))
anova(lmerTest::lmer(CopyNumber ~ type * Warm + (1|block) + (1|year), data = trait))

## GenomeSize
anova(lm(GenomeSize~group, trait))
otu_gsize_aov <- aov(lm(GenomeSize~group, trait))
TukeyHSD(otu_gsize_aov)
anova(lm(GenomeSize~type, trait))
t.test(subset(trait, type == "litter")$GenomeSize, subset(trait, type == "soil")$GenomeSize)

summary(lmerTest::lmer(GenomeSize ~ type * Warm + (1|block) + (1|year), data = trait))
anova(lmerTest::lmer(GenomeSize ~ type * Warm + (1|block) + (1|year), data = trait))

## GC
anova(lm(GC~group, trait))
otu_GC_aov <- aov(lm(GC~group, trait))
TukeyHSD(otu_GC_aov)
anova(lm(GC~type, trait))
t.test(subset(trait, type == "litter")$GC, subset(trait, type == "soil")$GC)

summary(lmerTest::lmer(GC ~ type * Warm + (1|block) + (1|year), data = trait))
anova(lmerTest::lmer(GC ~ type * Warm + (1|block) + (1|year), data = trait))

