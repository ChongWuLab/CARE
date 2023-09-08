library(tidyr)
library(ggplot2)
library(ggsci)

#setwd("/Users/cwu18/Dropbox/FSU_research/Undergoing/CARE/JASA-final/simulation-final/summaryRes/")
setwd("C:\\Users\\Evelyn\\Dropbox\\CARE\\simulation_final_2\\summaryRes\\")

files = list.files()
files = files[grepl("res_runtime_",files)]


out = NULL
for(j in 1:length(files)) {
    load(files[j])
    out = rbind(out,run_time2)
}

colMeans(out)
out = out[,c(-2,-3,-13)]
colnames(out) = c("CARE","cML","ContMix","MR-Lasso","IVW","MR-Egger","Weighted-Median","Weighted-Mode","RAPAS","MR-mix","MR-APSS")

out = out[,c(11,1,2,10,8,3,7,4,9,6,5)]
dat <- stack(as.data.frame(out))

ggplot(dat) +
  geom_boxplot(outlier.shape = NA,aes(x = ind, y = values)) + scale_y_continuous(limits = c(0,20)) + theme_bw() +
theme(axis.title = element_text(face="bold", size = rel(1.2)),
axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.text = element_text( size = rel(1.2)),
      axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
      legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) ) +  labs(
y = "Runtime (s)",
x = "Methods")+  scale_color_npg() + scale_fill_npg()

ggsave("runtime.jpeg")
