### subnetwork and multiple regression on matrices (MRM)
### Only codes for litterbag samples are shown, the others are similar.
library(readr)
library(ieggr)


############################
### 1. read in data
############################
# treatment
treat <- read.csv("data/treatment.csv")
rownames(treat) <- treat$sample
treat$group <- paste(treat$type, treat$Warm, sep = "_")
treat_L <- subset(treat, type == "litter")

GeoChip <- read_delim("data/GeoChip.txt")
GeoChip <- as.data.frame(GeoChip[,-1])
rownames(GeoChip) <- GeoChip$uniqueID

GeoChip_C <- subset(GeoChip, subcategory1 == "Carbon Degradation" | subcategory1 == "Carbon degradation")
unique(GeoChip_C$gene)
colnames(GeoChip_C)[11:ncol(GeoChip_C)] <- substring(colnames(GeoChip_C)[11:ncol(GeoChip_C)], 2, nchar(colnames(GeoChip_C)[11:ncol(GeoChip_C)]))
GeoChip_C_anno <- unique(GeoChip_C[, c(4,8,10)])




############################
### 2. response ratio
############################
x=1:5
y=9:13
x= x*0.0001
y = y * 0.0001
rr=respr(x,y,alternative = "two.sided")

### 2.1 litterbag
GeoChip_C_LW <- cbind(GeoChip_C[, 1:10], GeoChip_C[, colnames(GeoChip_C) %in% subset(treat_L, treat_L$Warm == "warming")$sample])
GeoChip_C_LC <- cbind(GeoChip_C[, 1:10], GeoChip_C[, colnames(GeoChip_C) %in% subset(treat_L, treat_L$Warm == "unwarming")$sample])

GeoChip_C_LW_2 <- GeoChip_C_LW[, 11:ncol(GeoChip_C_LW)]
GeoChip_C_LC_2 <- GeoChip_C_LC[, 11:ncol(GeoChip_C_LC)]

which(rownames(GeoChip_C_LW_2) != rownames(GeoChip_C_LC_2))

GeoChip_C_LW_3 <- GeoChip_C_LW_2[which(rowSums(GeoChip_C_LW_2, na.rm = T) > 0 | rowSums(GeoChip_C_LC_2, na.rm = T) > 0),]
GeoChip_C_LC_3 <- GeoChip_C_LC_2[which(rowSums(GeoChip_C_LW_2, na.rm = T) > 0 | rowSums(GeoChip_C_LC_2, na.rm = T) > 0),]

which(rownames(GeoChip_C_LW_3) != rownames(GeoChip_C_LC_3))

rr_L_WC <- as.data.frame(matrix(0, nrow = nrow(GeoChip_C_LW_3), ncol=15))

for (i in 1:nrow(rr_L_WC)) {
  rr_L_WC[i, ] = respr(GeoChip_C_LW_3[i, ], GeoChip_C_LC_3[i, ], alternative = "two.sided")
}
colnames(rr_L_WC) <- colnames(rr)
rownames(rr_L_WC) <- rownames(GeoChip_C_LW_3)

rr_L_WC$uniqueID <- rownames(rr_L_WC)
rr_L_WC <- merge(rr_L_WC, GeoChip_C, by = "uniqueID", all.y = F)

# rr_L_WC_sig <- subset(rr_L_WC, P < 0.05)
rr_L_WC_sig <- subset(rr_L_WC, P < 0.05 | RR == 'Inf' | RR == '-Inf')

# significant probes do response ratio again
rr_L_WC_sig2 <- rr_L_WC_sig
colnames(rr_L_WC_sig2)
rr_L_WC_sig2 <- aggregate(rr_L_WC_sig2[, 26:73], by=list(rr_L_WC_sig2$gene), FUN=sum, na.rm = T)
rownames(rr_L_WC_sig2) <- rr_L_WC_sig2$Group.1
rr_L_WC_sig2 <- rr_L_WC_sig2[, -1]

rr_L_WC_sig2_LW <- rr_L_WC_sig2[, colnames(rr_L_WC_sig2) %in% subset(treat_L, treat_L$Warm == "warming")$sample]
rr_L_WC_sig2_LC <- rr_L_WC_sig2[, colnames(rr_L_WC_sig2) %in% subset(treat_L, treat_L$Warm == "unwarming")$sample]
which(rownames(rr_L_WC_sig2_LW) != rownames(rr_L_WC_sig2_LC))

rr_L_WC_2 <- as.data.frame(matrix(0, nrow = nrow(rr_L_WC_sig2), ncol=15))
for (i in 1:nrow(rr_L_WC_2)) {
  rr_L_WC_2[i, ] = respr(rr_L_WC_sig2_LW[i, ], rr_L_WC_sig2_LC[i, ], alternative = "two.sided")
}
colnames(rr_L_WC_2) <- colnames(rr)
rownames(rr_L_WC_2) <- rownames(rr_L_WC_sig2)
