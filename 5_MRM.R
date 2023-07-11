### subnetwork and multiple regression on matrices (MRM)
### Only codes for bacterial network in warmed litterbag are shown, the others are similar.
setwd("D:/苏亦凡/My Libraries/litter bag/litterbag_GitHub/")
getwd()

library(reshape2)
library(vegan)
library(ieggr)
library(stringr)
library(igraph)
library(ecodist)


################################
### 1. read in data
################################
# treatment
treat <- read.csv("data/treatment.csv")
rownames(treat) <- treat$sample
treat$group <- paste(treat$type, treat$Warm, sep = "_")
treat_litter <- subset(treat, type == "litter")

# network
tranmatx = read.table("data/LitterWarming pairwise distance matrix.txt", header = T,row.names = 1)
tranmatx[tranmatx != 1] <- NA

otu_LW = read.table(paste0("data/otu_L_W.txt"), header = T, sep = "\t", row.names = 1)
colnames(otu_LW) <- str_sub(colnames(otu_LW), 3, str_count(colnames(otu_LW)))

property_16S <- read.table("data/subnetwork_property_16S.txt", sep = "\t", row.names = 1)
rownames(property_16S) <- property_16S$sample

# otu table
otu.file <- read.table("data/OTU_table_16S.txt", header = T, row.names = 1, check.names = F)
otu = as.data.frame(t(otu.file))
otu_litter <- subset(otu, rownames(otu) %in% treat_litter$sample)
otu_litter <- otu_litter[, colSums(otu_litter) != 0]

# environmental factors
env <- read.csv("data/env.csv", header = T)
env_litter <- env
rownames(env_litter) <- env_litter$sample
env_litter <- env_litter[,c(1:5, 7, 14:24, 27:33)]



########################################################
### 2. calculate topological properties for subnetwork
########################################################
#transition matrix
adj_LW = tranmatx
otu_LW = t(otu_LW)

#match names#
match.rowname = match.name(both.list = list(adj_LW=adj_LW),cn.list = list(otu_LW=otu_LW))
otu_LW = match.rowname$otu_LW
dim(otu_LW)
adj_LW=match.rowname$adj_LW
dim(adj_LW)

#transfer adj_LW into matrix for igraph
m_LW = as.matrix(adj_LW)

#unique otu names appeared in each sample
otu_LW = t(otu_LW)
id_LW = sapply(1:ncol(otu_LW), function(i){rownames(otu_LW[which(otu_LW[,i]!=0),])})

#create global network #
my_graph_LW = graph_from_adjacency_matrix(m_LW, mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL)

#subnetwork property#
sample.out_LW = matrix(0,ncol(otu_LW),8)
for(k in 1:ncol(otu_LW)){
  g <- subgraph(my_graph_LW, v=id_LW[[k]])# this generate a subnetwork only containing genes/OTUs in a certain sample#

  sample.out_LW[k,1] = gorder(g) 
  sample.out_LW[k,2] = ecount(g)
  sample.out_LW[k,3] = mean(centr_degree(g)$res)
  sample.out_LW[k,4] = avgCC = transitivity(g, type="average", isolates="zero")
  sample.out_LW[k,5] = mean_distance(g, directed=F, unconnected=T)
  
  gd = cluster_fast_greedy(g)
  sample.out_LW[k,6] = length(gd)
  sample.out_LW[k,7] = modularity(gd)
  sample.out_LW[k,8] = graph.density(g)
}
colnames(sample.out_LW) = c("node number", "link number", "avgK", "avgCC", 
                             "geodesic distance", "number of modules", "Modularity",
                             "density")
sample.out_LW <- as.data.frame(sample.out_LW)
sample.out_LW$sample = colnames(otu_LW)

property_LW <- sample.out_LW
rownames(property_LW) <- property_LW$sample




#################
### 3. MRM
#################
match.rowname=match.name(rn.list = list(env_litter = env_litter, otu_litter = otu_litter, property_16S = property_16S))
property_16S_L = match.rowname$property_16S
env_litter = match.rowname$env_litter
otu_litter = match.rowname$otu_litter

dim(env_litter);dim(property_16S_L)

# community
dis_16S_L <- vegdist(otu_litter, method = "bray")

# biotic
node.number_16S_L <- vegdist(scale(property_16S_L$node.number),"euclid")
link.number_16S_L <- vegdist(scale(property_16S_L$link.number),"euclid")
avgK_16S_L <- vegdist(scale(property_16S_L$avgK),"euclid")
avgCC_16S_L <- vegdist(scale(property_16S_L$avgCC),"euclid")
geodesic.distance_16S_L <- vegdist(scale(property_16S_L$geodesic.distance),"euclid")
number.of.modules_16S_L <- vegdist(scale(property_16S_L$number.of.modules),"euclid")
Modularity_16S_L <- vegdist(scale(property_16S_L$Modularity),"euclid")
density_16S_L <- vegdist(scale(property_16S_L$density),"euclid")

# abiotic
pH_L <- vegdist(scale(env_litter$pH),"euclidean")
FlTotl_L <- vegdist(scale(env_litter$FlTotl),"euclidean")
FlC3_L <- vegdist(scale(env_litter$FlC3),"euclidean")
FlC4_L <- vegdist(scale(env_litter$FlC4),"euclidean")
plant.richness_L <- vegdist(scale(env_litter$plant.richness),"euclidean")
moisture_samplingmonth_L <- vegdist(scale(env_litter$moisture_samplingmonth),"euclidean")
annual_moisture_L <- vegdist(scale(env_litter$annual_moisture),"euclidean")
temperature_annual_L <- vegdist(scale(env_litter$temperature_annual),"euclidean")
ER_annualmean_L <- vegdist(scale(env_litter$ER_annualmean),"euclidean")
GPP_annualmean_L <- vegdist(scale(env_litter$GPP_annualmean),"euclidean")
NEE_annualmean_L <- vegdist(scale(env_litter$NEE_annualmean),"euclidean")
Autotrophic_L <- vegdist(scale(env_litter$Autotrophic),"euclidean")
Heterotrophic_L <- vegdist(scale(env_litter$Heterotrophic),"euclidean")
total_soil_respiration_L <- vegdist(scale(env_litter$total_soil_respiration),"euclidean")


# select variables and compare different model
# only the final models are shown here
# MRM of abiotic factors
MRM_abiotic_16S_L = MRM(dis_16S_L ~ moisture_samplingmonth_L + FlC4_L +  temperature_annual_L, 
                        nperm=1000)
MRM_abiotic_16S_L

# MRM of biotic factors
MRM_biotic_16S_L = MRM(dis_16S_L ~ link.number_16S_L, nperm=1000)
MRM_biotic_16S_L

# MRM of all factors #
MRM_all_16S_L = MRM(dis_16S_L ~ moisture_samplingmonth_L + FlC4_L +  temperature_annual_L + link.number_16S_L,nperm=1000)
MRM_all_16S_L

