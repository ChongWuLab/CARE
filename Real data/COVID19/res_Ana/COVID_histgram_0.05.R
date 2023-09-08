library(ggplot2)
library(dplyr)
library(tidyr)
# prepare the data
#setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/CARE/RealData/New/Res_Ana/")
setwd("/rsrch5/scratch/biostatistics/wzhang24/CARE/covid19/res_Ana/")
load("COVID19_summarized.RData")

# All the used exposures; others are duplicated...

#used.id = c("ebi-a-GCST006099","ieu-b-24","ieu-b-73","ieu-b-2","ebi-a-GCST90000025","ebi-a-GCST008403","ieu-a-44","ebi-a-GCST004634","ebi-a-GCST004631","ieu-b-41","ieu-a-1083","ieu-a-1096","ieu-a-1102","ebi-a-GCST003837","ieu-a-836","ebi-a-GCST005195","ieu-a-12","ieu-b-39","ebi-a-GCST010723","ieu-a-996","ebi-a-GCST006464","ebi-a-GCST004606","ebi-a-GCST004617","ebi-a-GCST004600","ebi-a-GCST006944","ieu-b-109","ebi-a-GCST009541","ebi-a-GCST006979","ieu-a-89","ebi-a-GCST004604","ieu-a-31","ebi-a-GCST006908","ebi-a-GCST007090","ieu-b-110","ieu-a-966","ebi-a-GCST004625","ebi-a-GCST004609","ieu-b-18","ebi-a-GCST004626","ieu-a-798","ebi-a-GCST006940","ebi-a-GCST004629","ebi-a-GCST004623","ebi-a-GCST004633","ebi-a-GCST007092","ebi-a-GCST007091","ieu-b-7","ebi-a-GCST007430","ebi-a-GCST004616","ebi-a-GCST006085","ebi-a-GCST004601","ebi-a-GCST006804","ebi-a-GCST004619","ieu-a-833","ieu-a-22","ebi-a-GCST90000616","ebi-a-GCST003839","ebi-a-GCST009760","ieu-a-967","ebi-a-GCST006100","ebi-a-GCST006906","ebi-a-GCST004621","ebi-a-GCST004613","ieu-b-38","ieu-a-301","ebi-a-GCST006867","ieu-a-970","ebi-a-GCST006586","ieu-a-60","ieu-a-72","ebi-a-GCST004610","ebi-a-GCST004618","ebi-a-GCST003769","ebi-a-GCST005902","ieu-a-1000","ebi-a-GCST005038","ebi-a-GCST004988","ieu-a-7","ieu-a-85","ieu-a-86","ieu-a-87","ebi-a-GCST009722","ebi-a-GCST004614","ebi-a-GCST004608","ebi-a-GCST004615","ebi-a-GCST004611","ebi-a-GCST004628","ebi-a-GCST90002412","ebi-a-GCST007431","ebi-a-GCST007429","ebi-a-GCST004627","ebi-a-GCST004632","ebi-a-GCST004630","ebi-a-GCST004605","ebi-a-GCST004602","ebi-a-GCST004599","ieu-a-93","ieu-a-975","ebi-a-GCST006097","ebi-a-GCST006098","ieu-b-111","ebi-a-GCST006947","ebi-a-GCST006945","ebi-a-GCST006942","ebi-a-GCST006951","ebi-a-GCST006943","ebi-a-GCST006948","ebi-a-GCST006952","ebi-a-GCST006950","ebi-a-GCST006941","ebi-a-GCST006478","ebi-a-GCST006946","ebi-a-GCST90000027","ebi-a-GCST005920","ebi-a-GCST005921","ebi-a-GCST005923","ebi-a-GCST006696","ebi-a-GCST006697","ebi-a-GCST006368","ebi-a-GCST006698","ieu-a-90","ieu-a-91","ebi-a-GCST006475","ieu-a-965")

#CARE = CARE.res[CARE.res[,"exposure.id"] %in% used.id,]
CARE = CARE.res
#colSums(CARE.res<0.05/140,na.rm =T)
#CARE.res[CARE.res[,"BIC.p"]<0.05/140,]
#colSums(CARE<0.01,na.rm=T)


CARE = CARE[,c("Name","BIC.tilde_theta","BIC.se","BIC.p","exposure.id","outcome.id","exposure.name","outcome.name")]

colnames(CARE) = c("Trait","CARE_beta","CARE_se","CARE_p","exposure.id","outcome.id","exposure.name","outcome.name")

#sum(CARE$CARE_p<0.05/137)


MR.res2 = MR.res[,c("DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","IVW_beta","IVW_se","IVW_p","RAPS_beta","RAPS_se","RAPS_p","MRmix_beta","Weighted_median_beta","Weighted_median_p","MRAPSS_beta","MRAPSS_se","MRAPSS_p","Trait")]

dat = left_join(CARE,MR.res2,by = "Trait")

dat = left_join(CARE,MR.res,by = "Trait")

dat$outcome.name = dat$outcome.id #outcome id is empty;



cutoff = 0.05#0.05/140
colSums(dat<0.05/140,na.rm=T)

colSums(dat<0.01,na.rm=T)

#,"IVW_fixed_p"
dat1 = dat[dat[,"outcome.id"]=="B2",]

dat1[dat1$CARE_p<cutoff,]
B2 = colSums(dat1< cutoff,na.rm=T)

B2 = B2[c("CARE_p","DP_MA_BIC_p","ConMix_p","MRLasso_p","IVW_p","Egger_p","Weighted_median_p","Weighted_mode_p","RAPS_p","MRmix_p","MRAPSS_p")]


#dat2 = dat[dat[,"outcome.id"]=="C2",]

#C2 =  colSums(dat2< cutoff,na.rm=T)
#C2 = C2[c("CARE_p","DP_MA_BIC_p","ConMix_p","MRLasso_p","IVW_p","Egger_p","Weighted_median_p","Weighted_mode_p","RAPS_p","MRmix_p","MRAPSS_p")]

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

pdf(paste("COVID19_barplot_0.05.pdf"),
    paper = "special",
    height = 8,width = 11)
fp1
dev.off()

jpeg("COVID19_barplot_0.05.jpeg",width=9,height=6,units="in",res=300)
fp1
dev.off()
