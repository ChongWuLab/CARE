library(ggplot2)
library(dplyr)
library(tidyr)
# prepare the data

setwd("/rsrch5/scratch/biostatistics/wzhang24/CARE/covid19/res_Ana/")
load("COVID19_summarized.RData")


CARE = CARE.res

CARE = CARE[,c("Name","BIC.tilde_theta","BIC.se","BIC.p","exposure.id","outcome.id","exposure.name","outcome.name")]

colnames(CARE) = c("Trait","CARE_beta","CARE_se","CARE_p","exposure.id","outcome.id","exposure.name","outcome.name")


MR.res2 = MR.res[,c("DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","IVW_beta","IVW_se","IVW_p","RAPS_beta","RAPS_se","RAPS_p","MRmix_beta","Weighted_median_beta","Weighted_median_p","MRAPSS_beta","MRAPSS_se","MRAPSS_p","Trait")]

dat = left_join(CARE,MR.res2,by = "Trait")

dat = left_join(CARE,MR.res,by = "Trait")

dat$outcome.name = dat$outcome.id #outcome id is empty;



cutoff = 0.05/45


#,"IVW_fixed_p"
dat1 = dat[dat[,"outcome.id"]=="B2",]

dat1[dat1$CARE_p<cutoff,]
B2 = colSums(dat1< cutoff,na.rm=T)

B2 = B2[c("CARE_p","DP_MA_BIC_p","ConMix_p","MRLasso_p","IVW_p","Egger_p","Weighted_median_p","Weighted_mode_p","RAPS_p","MRmix_p","MRAPSS_p")]

#dat = rbind(B2,C2)
#dat = t(dat)
dat = B2
dat = as.data.frame(dat)
colnames(dat) = c("B2: Hospitalized COVID-19 vs population")
rownames(dat) = c("CARE","cML-DP","ContMix","MR-Lasso","IVW","MR-Egger","Median", "Mode","RAPS","MR-mix","MR-APSS")

dat$methods = rownames(dat)

final_dat = dat %>% gather("Diseases","Value",-methods)


library("ggsci")
pal_npg("nrc")(3)

library(scales)

fp1 = ggplot(data=final_dat, aes(x=reorder(methods,-Value), y=Value)) +
  geom_bar(stat="identity", fill = "#74AED4", show.legend = FALSE) +
  xlab("Methods") +
  ylab("# of significant risk factors") +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.grid.major.y = element_line(color = "#00A087FF", size = 0.5, linetype = 2)
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(breaks= pretty_breaks())

#pdf(paste("COVID19_barplot_bonferroni.pdf"),
#    paper = "special",
#    height = 8,width = 11)
#fp1
#dev.off()

jpeg("COVID19_barplot_bonferroni.jpeg",width=9,height=6,units="in",res=300)
fp1
dev.off()

