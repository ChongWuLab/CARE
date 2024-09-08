library(ggplot2)
library(dplyr)
library(tidyr)
# prepare the data
setwd("CARE/covid19/res_Ana/")
load("COVID19_summarized.RData")

# All the used exposures; others are duplicated...
#used.id = c("ebi-a-GCST006099","ieu-b-24","ieu-b-73","ieu-b-2","ebi-a-GCST90000025","ebi-a-GCST008403","ieu-a-44","ebi-a-GCST004634","ebi-a-GCST004631","ieu-b-41","ieu-a-1083","ieu-a-1096","ieu-a-1102","ebi-a-GCST003837","ieu-a-836","ebi-a-GCST005195","ieu-a-12","ieu-b-39","ebi-a-GCST010723","ieu-a-996","ebi-a-GCST006464","ebi-a-GCST004606","ebi-a-GCST004617","ebi-a-GCST004600","ebi-a-GCST006944","ieu-b-109","ebi-a-GCST009541","ebi-a-GCST006979","ieu-a-89","ebi-a-GCST004604","ieu-a-31","ebi-a-GCST006908","ebi-a-GCST007090","ieu-b-110","ieu-a-966","ebi-a-GCST004625","ebi-a-GCST004609","ieu-b-18","ebi-a-GCST004626","ieu-a-798","ebi-a-GCST006940","ebi-a-GCST004629","ebi-a-GCST004623","ebi-a-GCST004633","ebi-a-GCST007092","ebi-a-GCST007091","ieu-b-7","ebi-a-GCST007430","ebi-a-GCST004616","ebi-a-GCST006085","ebi-a-GCST004601","ebi-a-GCST006804","ebi-a-GCST004619","ieu-a-833","ieu-a-22","ebi-a-GCST90000616","ebi-a-GCST003839","ebi-a-GCST009760","ieu-a-967","ebi-a-GCST006100","ebi-a-GCST006906","ebi-a-GCST004621","ebi-a-GCST004613","ieu-b-38","ieu-a-301","ebi-a-GCST006867","ieu-a-970","ebi-a-GCST006586","ieu-a-60","ieu-a-72","ebi-a-GCST004610","ebi-a-GCST004618","ebi-a-GCST003769","ebi-a-GCST005902","ieu-a-1000","ebi-a-GCST005038","ebi-a-GCST004988","ieu-a-7","ieu-a-85","ieu-a-86","ieu-a-87","ebi-a-GCST009722","ebi-a-GCST004614","ebi-a-GCST004608","ebi-a-GCST004615","ebi-a-GCST004611","ebi-a-GCST004628","ebi-a-GCST90002412","ebi-a-GCST007431","ebi-a-GCST007429","ebi-a-GCST004627","ebi-a-GCST004632","ebi-a-GCST004630","ebi-a-GCST004605","ebi-a-GCST004602","ebi-a-GCST004599","ieu-a-93","ieu-a-975","ebi-a-GCST006097","ebi-a-GCST006098","ieu-b-111","ebi-a-GCST006947","ebi-a-GCST006945","ebi-a-GCST006942","ebi-a-GCST006951","ebi-a-GCST006943","ebi-a-GCST006948","ebi-a-GCST006952","ebi-a-GCST006950","ebi-a-GCST006941","ebi-a-GCST006478","ebi-a-GCST006946","ebi-a-GCST90000027","ebi-a-GCST005920","ebi-a-GCST005921","ebi-a-GCST005923","ebi-a-GCST006696","ebi-a-GCST006697","ebi-a-GCST006368","ebi-a-GCST006698","ieu-a-90","ieu-a-91","ebi-a-GCST006475","ieu-a-965") #,"ebi-a-GCST006368"; remove this, this is duplicated BMI


#CARE = CARE.res[CARE.res[,"exposure.id"] %in% used.id,]
CARE = CARE.res
CARE = CARE[,c("Name","BIC.tilde_theta","BIC.se","BIC.p","exposure.id","outcome.id","exposure.name","outcome.name")]

colnames(CARE) = c("Trait","CARE_beta","CARE_se","CARE_p","exposure.id","outcome.id","exposure.name","outcome.name")

#sum(CARE$CARE_p<0.05/1)
#MR.res2 = MR.res[,c("DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","IVW_beta","IVW_se","IVW_p","RAPS_beta","RAPS_se","RAPS_p","MRmix_beta","Weighted_median_beta","Weighted_median_p","Trait")]
MR.res2 = MR.res[,c("DP_MA_BIC_beta","DP_MA_BIC_se","DP_MA_BIC_p","IVW_beta","IVW_se","IVW_p","MRAPSS_beta","MRAPSS_se","MRAPSS_p","Trait")]

dat = left_join(CARE,MR.res2,by = "Trait")
dat$outcome.name = dat$outcome.id #outcome id is empty;


dat[is.na(dat)] = 1
#cutoff = 0.001
#indx = dat[,"CARE_p"]< cutoff | dat[,"DP_MA_BIC_p"]< cutoff |dat[,"IVW_p"]< cutoff |dat[,"RAPS_p"]< cutoff |dat[,"Weighted_median_p"]< cutoff
#tmp = dat[indx,]
tmp = dat[dat[,"outcome.id"] == "B2",]

indx = p.adjust(tmp$CARE_p)<0.05 | p.adjust(tmp$DP_MA_BIC_p)<0.05 | p.adjust(tmp$IVW_p)<0.05 | p.adjust(tmp$MRAPSS_p)<0.05
tmp = tmp[indx,]

sig.exp = tmp[,"exposure.id"]

final_plot_data = dat[dat[,"exposure.id"] %in% sig.exp,]
final_plot_data = final_plot_data[,c("DP_MA_BIC_beta","DP_MA_BIC_p","IVW_beta","IVW_p","MRAPSS_beta","MRAPSS_p","CARE_beta","CARE_p","exposure.name","outcome.name")]

final_plot_data[final_plot_data[,"outcome.name"] == "B2","outcome.name"] = "B2: Hospitalized covid vs. population"
#final_plot_data[final_plot_data[,"outcome.name"] == "C2_all_eur_V6","outcome.name"] = "C2: Covid vs. population"

#final_plot_data[final_plot_data[,"exposure.name"] == "Family history of Alzheimer's disease","outcome.name"] =

final_plot_data = final_plot_data[final_plot_data['outcome.name']=="B2: Hospitalized covid vs. population",]

final_dat = final_plot_data %>% gather("Methods","Value",-exposure.name,-outcome.name)
final_dat$Method = gsub("_beta","",final_dat$Methods)
final_dat$Method = gsub("_p","",final_dat$Method)
final_dat$Methods = gsub(".*_","",final_dat$Methods)

final_dat[final_dat$Method=="DP_MA_BIC","Method"] = "cML-DP"
final_dat[final_dat$Method=="MRAPSS","Method"] = "MR-APSS"

final_dat = final_dat %>% spread(Methods,Value)

final_dat = final_dat[,c(3,5,4,1,2)]
colnames(final_dat) = c("Method","pval","est","trait1","trait2")
#Method          pval           est  trait1 trait2

#final_dat$trait1 = factor(final_dat$trait1, levels=c('Body mass index','Waist circumference','Obesity class 1','Experiencing mood swings','Depressed affect'))#'Overweight',,'multiple sclerosis'
final_dat$trait1 = factor(final_dat$trait1)


# start plot --------------------------------------------------------------
#finish the plot


### plot data
all_plot_data = final_dat
all_plot_data$Method = factor(all_plot_data$Method,
unique(all_plot_data$Method))

pointsize0 = -log10(all_plot_data$pval)

pointsize = pointsize0
pointsize[pointsize0 < -log10(0.05)] = 1
pointsize[(pointsize0 > -log10(0.05)) & (pointsize0 < -log10(0.05/50))] = 2
pointsize[pointsize0 > -log10(0.05/50)] = 3
pointsize = as.factor(pointsize)
pointsize = factor(pointsize,
labels = c("p > 0.05",
"0.001 < p < 0.05",
"p < 0.001"))
all_plot_data = cbind(all_plot_data,pointsize)

pointshape0 = all_plot_data$est > 0
pointshape = pointshape0
pointshape[pointshape0] = "+"
pointshape[!pointshape0] = "-"
all_plot_data = cbind(all_plot_data,pointshape)

pointfill = rep(0,nrow(all_plot_data))
pointfill[pointsize == "p < 0.001"] = 0
pointfill[pointsize != "p < 0.001"] = 1
#pointfill = as.factor(pointfill)
all_plot_data = cbind(all_plot_data,pointfill)

x_axis = as.numeric(all_plot_data$Method)
all_plot_data = cbind(all_plot_data,x_axis)


# using the color suggested
library("ggsci")
pal_npg("nrc")(5)

### plot
fp1 = ggplot(all_plot_data, aes(x=Method, y=1)) +
geom_point(aes(color = Method,
size = pointsize,
shape = pointshape,
fill = factor(ifelse(pointfill, "No", Method))
)) +
facet_grid(trait1~trait2,switch = "y") +
theme(strip.background = element_blank(),
legend.title = element_blank(),
legend.position = "bottom",
legend.direction="vertical",
axis.text.x = element_blank(),
axis.text.y = element_blank(),
axis.ticks = element_blank(),
panel.spacing.x=unit(0, "lines"),
panel.spacing.y=unit(0, "lines"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.background = element_blank(),
panel.border = element_rect(colour = "black",
fill = NA,
size = 0.5),strip.text.y.left = element_text(angle = 0),strip.text.y = element_text(size = 12),strip.text.x = element_text(size = 12)) +
labs(y = "", x = "") + scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","white"), na.value=NA, guide="none") +
scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")) + scale_shape_manual(values=c(25,24)) +
scale_size_manual(values = c(1,2.5,5)) +
guides(colour = guide_legend(override.aes = list(shape = 17,size = 2.5)),
size = guide_legend(override.aes = list(shape = c(2,2,17))))


pdf(paste("COVID19.pdf"),
paper = "special",
height = 8,width = 8)
fp1
dev.off()

jpeg("COVID19_point.jpeg",height=8,width=8,units="in",res=300)
fp1
dev.off()


strip.text.y = element_text(angle = 0)

+
scale_fill_manual(values=c("#6A3D9A","#FDBF6F","#1F78B4", "#E31A1C",
"white", "#FB9A99","#B2DF8A",
"#FF7F00", "#CAB2D6",
"#A6CEE3"),
na.value=NA, guide="none",
breaks=c("Causal",
"Correlated",
"Unrelated",
"Non-causal")) +
scale_color_manual(values=c("#6A3D9A","#FDBF6F","#1F78B4", "#E31A1C",
"#FB9A99","#B2DF8A",
"#FDBF6F", "#FF7F00", "#CAB2D6",
"#A6CEE3")) 
