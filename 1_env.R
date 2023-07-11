### linear mixed model for environmental factors
library(reshape2)


##########################
### 1. read in data
##########################
env <- read.csv("data/env.csv", header = T)
loss <- read.table("data/mass_loss.txt", header = T, sep = ";")



###################################
### 2. warming effect on mass loss
###################################
lm1 <- lmerTest::lmer(loss ~ warm + (1|block) + (1|year), data=loss)
summary(lm1)
anova(lm1)

summary(lmerTest::lmer(loss ~ warm + (1|block), data=subset(loss, year == 2016)))
summary(lmerTest::lmer(loss ~ warm + (1|block), data=subset(loss, year == 2017)))




#####################################################
### 3. linear mixed model for environmental factors
#####################################################
env_name <- colnames(env)[c(14:24, 27:33)]
env2 <- melt(env, measure.vars = env_name)

env_LMM <- as.data.frame(matrix(nrow = length(env_name), ncol=6))
colnames(env_LMM) <- c("env", "estimate", "std.error", "df", "t_value", "Pr(>|t|)")

for (i in 1:length(env_name)) {
  subenv <- subset(env2, variable == env_name[i])
  lm <- lmerTest::lmer(value ~ Warm + (1|block) + (1|year), data = subenv)
  lm_summ <- summary(lm)
  env_LMM[i, "env"] <- env_name[i]
  env_LMM[i, "estimate"] <- lm_summ$coefficients[2,1]
  env_LMM[i, "std.error"] <- lm_summ$coefficients[2,2]
  env_LMM[i, "df"] <- lm_summ$coefficients[2,3]
  env_LMM[i, "t_value"] <- lm_summ$coefficients[2,4]
  env_LMM[i, "Pr(>|t|)"] <- lm_summ$coefficients[2,5]
}

env_LMM$padj <- p.adjust(env_LMM$`Pr(>|t|)`, method = 'BH') 




