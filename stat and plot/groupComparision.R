#########################
library(readxl)
library("diffcyt")
library(lme4)
library(lmerTest)
library("multcomp")

types = c("Fig2b", "Fig1F NCa_DN", "Fig1F DN_T", "Fig 2A", "Fig 1E NCa vs DN", "Fig 1E DN vs veHCC", "Fig S6C", "Fig S6F", "Fig 6G oncogene", "Fig S6G TSG", "DN VS veHCC-CNA", "DN VS veHCC-CFG")
p_vall1 = rep(NA, length(types))
names(p_vall1) = types
p_vall2 = p_vall1 
for(s in types ){
  data = as.matrix(data.frame(read_excel("SNV-TERT expression-LMM.xlsx", sheet=s, col_names=F)))
  group = sub("\\d", "", sub(".*_", "", sub("_R.*", "", data[,1])))
  print(table(group))
  data1 = data.frame( y=as.numeric(data[,2]), group=group, sample=data[,1], patient=sub("_.*", "", data[,1]) )
  fit <- lmer(y ~ group + (1 | patient), data = data1)
  # t 检验（Satterthwaite/Kenward-Roger 法）
  p_vall2[s] <- summary(fit)$coefficients[2, "Pr(>|t|)"]
  # Wald 检验（基于渐近正态分布）
  contrast =  matrix(c(0, 1), nrow = 1, ncol = 2,  dimnames = list(NULL, c("(Intercept)", "groupDN")) )
  test <- glht(fit, contrast)
  p_vall1[s] <- summary(test)$test$pvalues
}
write.csv( cbind( p_vall1,  p_vall2 ), "P.csv")

