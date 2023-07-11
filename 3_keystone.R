### Identify keystones
### Only codes for 16S are shown, the codes for ITS are similar.

library(edgeR)
library(ieggr)
library(vegan)
library(ALDEx2)
library(dplyr)

#####################
### 1. read in data
#####################
otu.file <- read.table("data/OTU_table_16S.txt", header = T, row.names = 1, check.names = F)
ITS.file <- read.table("data/OTU_table_ITS.txt", header = T, row.names = 1, check.names = F)
treat <- read.csv("data/treatment.csv")
rownames(treat) <- treat$sample

otu = as.matrix(otu.file)
ITS = as.matrix(ITS.file)



#####################
### 2. edgeR
#####################
scp <- match.name(rn.list = list(treat = treat), cn.list = list(otu = otu))
otu2 <- scp$otu
treat <- scp$treat
group1 <- treat$type

dgelist1 <- DGEList(counts = otu2, group = treat$type)

# 2.1 filter low count
# a.present in at least 1/4 samples
# b.relative abundance >= 0.1% in at least 1/8 samples
otu2 <- as.data.frame(otu2)
otu2_ra <- decostand(otu2, "total", MARGIN = 2)
colSums(otu2_ra)
otu2_ra[otu2_ra >= 0.001] <- 2
otu2_ra[otu2_ra > 0 & otu2_ra < 0.001] <- 1

otu2_ra <- otu2_ra[rowSums(otu2_ra == 0) <= 12, ]
otu2_ra <- otu2_ra[rowSums(otu2_ra == 2) >= 6, ]

otu2 <- otu2[rownames(otu2_ra),]
scp <- match.name(rn.list = list(treat = treat), cn.list = list(otu2 = otu2))

dgelist1 <- DGEList(counts = otu2, group = treat$type)

# 2.2 standardize
dgelist1_norm <- calcNormFactors(dgelist1, method = 'TMM')
plotMDS(dgelist1_norm)

design1 <- model.matrix(~group1)
design1[design1[,2]==0, 2] <- 2
design1[design1[,2]==1, 2] <- 0
design1[design1[,2]==2, 2] <- 1

dge1 <- estimateDisp(dgelist1_norm, design1, robust = TRUE)
plotBCV(dge1)

# 2.3 quasi-likelihood negative binomial generalized log-linear model 拟合
fit1 <- glmQLFit(dge1, design1, robust = TRUE) #拟合模型
lrt1 <- glmQLFTest(fit1)    #统计检验

topTags(lrt1)

output1 <- as.data.frame(topTags(lrt1, n = nrow(dgelist1$counts)))
output1$ID <- rownames(output1)
output1 <- merge(output1, AveLogCPM1, by = "ID")

# 2.4 output significant results
edgeR_16S <- output1
edgeR_16S[which(edgeR_16S$FDR < 0.01 & edgeR_16S$logFC <= -1),'sig'] <- 'depleted'
edgeR_16S[which(edgeR_16S$FDR < 0.01 & edgeR_16S$logFC >= 1),'sig'] <- 'enriched'
edgeR_16S[which(edgeR_16S$FDR >= 0.01 | abs(edgeR_16S$logFC) < 1),'sig'] <- 'no diff'



#####################
### 3. ALDEx2
#####################
x_16S <- aldex.clr(otu2, group1, mc.samples=128, verbose=TRUE)
x.tt_16S <- aldex.ttest(x_16S, group1, paired.test=FALSE)
head(x.tt_16S)
x.effect_16S <- aldex.effect(x_16S, verbose=FALSE)
x.all_16S <- data.frame(x.tt_16S, x.effect_16S)
x.all_16S <- aldex(otu2, group1, mc.samples=128, test="t", effect=TRUE,
                   include.sample.summary=FALSE, denom="all", verbose=FALSE)
aldex.plot(x.all_16S, type="MA", test="welch", xlab="Log-ratio abundance",ylab="Difference")

sig_both_16S <- subset(x.all_16S, x.all_16S$we.eBH < 0.05 & x.all_16S$wi.eBH < 0.05)
edgeR_16S_sig <- subset(edgeR_16S, sig != "no diff")
length(which(rownames(sig_both_16S) %in% edgeR_16S_sig$ID))



#####################
### 4. overlap
#####################
sig_both_16S2 <- sig_both_16S
sig_both_16S2$ID <- rownames(sig_both_16S2)
sig_both_16S2[which(sig_both_16S2$effect > 0),'sig'] <- 'depleted'
sig_both_16S2[which(sig_both_16S2$effect < 0),'sig'] <- 'enriched'
sig_both_16S2 <- sig_both_16S2[, 12:13]

sig_16S <- as.data.frame(matrix(edgeR_16S$ID))
colnames(sig_16S) <- "ID"
for (i in sig_16S$ID) {
  if(!(i %in% sig_both_16S2$ID)){
    sig_16S[sig_16S$ID==i,'sig'] <- 'no diff'
  }
  else if(edgeR_16S[edgeR_16S$ID==i, "sig"] == "enriched" & sig_both_16S2[sig_both_16S2$ID==i, "sig"] == "enriched"){
    sig_16S[sig_16S$ID==i,'sig'] <- 'enriched'
  }
  else if(edgeR_16S[edgeR_16S$ID==i, "sig"] == "depleted" & sig_both_16S2[sig_both_16S2$ID==i, "sig"] == "depleted"){
    sig_16S[sig_16S$ID==i,'sig'] <- 'depleted'
  }
  else if(edgeR_16S[edgeR_16S$ID==i, "sig"] == "no diff"){
    sig_16S[sig_16S$ID==i,'sig'] <- 'no diff'
  }
}
sig_16S %>% group_by(sig) %>% summarise(n())

