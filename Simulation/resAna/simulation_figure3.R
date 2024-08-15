library(tidyr)
library(ggplot2)
library(ggsci)

#setwd("/Users/cwu18/Dropbox/FSU_research/Undergoing/CARE/JASA-final/simulation-final/summaryRes/")
setwd("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_final_2/summaryRes/")
files = list.files(getwd())
#ind = 'simRes1'
#files = files[grepl(ind,files)]
#ind = 'prop_invalid0.3'
#files = files[grepl(ind,files)]
draw.figure <- function(figureTitle,saveTitle,MSE_trunc =0.003) {
    #parameter:
    #MSE_trunc = 0.003
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        #file = "CD2_simRes11theta-0.2prop_invalid0.3.RData"
        load(file)
        final.out
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,1],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    df = df[,c("MRcML","MRcML_DP","IVW","Weighted_Median","Weighted_Mode","CARE_Efron_BIC","MRAPSS","theta")]
    
    colnames(df) = c("cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","CARE","MRAPSS","theta")
    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","MRAPSS"))
    # Power # Bias
    
    #
    p1 = ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) + geom_line(aes(linetype = Method)) +
    geom_point(aes(shape =Method), size = 2)+geom_hline(yintercept = 0.05) +   scale_y_continuous(breaks = sort(c(seq(0,1, length.out=5), 0.05))) + scale_shape_manual(values =c(2,3,4,5,6,7,8)) +  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash','solid')) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
    labs(title=paste0(figureTitle," (A)") ,
         y = "Power",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + guides(col =
         guide_legend(nrow = 1),shape =guide_legend(nrow = 1) )  + scale_color_npg() + scale_fill_npg()
         
    
    
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(abs(final.out[,"Bias"]),tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    df = df[,c("MRcML","MRcML_DP","IVW","Weighted_Median","Weighted_Mode","CARE_Efron_BIC","MRAPSS","theta")]
    
    df[,1:7] = abs(df[,1:7])#
    #df = df[df[,"theta"]!=0,]
    colnames(df) = c("cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","CARE","MRAPSS","theta")
    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    df.plot[df.plot[,3]>1,3] = 1

    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","MRAPSS"))
    # Power # Bias
    
    
    #stat_smooth(aes(fill=Method, linetype=Method), size=0.5, se = FALSE, level = 0.95) +
    p2 = ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) + geom_line(aes(linetype = Method)) +  geom_point(aes(shape =Method), size = 2) + scale_shape_manual(values =c(2,3,4,5,6,7,8)) +  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash',"solid")) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
          labs(title=paste0(figureTitle," (B)") ,
         y = "Absolute Bias",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_npg() + scale_fill_npg()

    
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,"MSE"],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    df = data.frame(out)
    
    df = df[,c("MRcML","MRcML_DP","IVW","Weighted_Median","Weighted_Mode","CARE_Efron_BIC","MRAPSS","theta")]
    
    colnames(df) = c("cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","CARE","MRAPSS","theta")
    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    
    df.plot[df.plot[,"Power"]>MSE_trunc,"Power"] = MSE_trunc

    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","MRAPSS"))
    
    #
    p3 = ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) +geom_line(aes(linetype = Method)) +  geom_point(aes(shape =Method), size = 1.5, alpha=0.85)  + scale_shape_manual(values =c(2,3,4,5,6,7,8)) +   scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash','solid')) + ylim(0,MSE_trunc) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
    axis.text = element_text( size = rel(1.2)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
    legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
    labs(title=paste0(figureTitle," (C)") ,
    y = "Mean Squared Error",
    x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_npg() + scale_fill_npg()
    
    
    
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,"Coverage"],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    df = data.frame(out)
    
    df = df[,c("MRcML","MRcML_DP","IVW","Weighted_Median","Weighted_Mode","CARE_Efron_BIC","MRAPSS","theta")]
    
    colnames(df) = c("cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","CARE","MRAPSS","theta")
    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","cML","cML-DP","IVW","Weighted-Median","Weighted-Mode","MRAPSS"))
    
    df.plot[df.plot[,3]<=0.6,3] = 0.6

    #
    p4 = ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) +geom_line(aes(linetype = Method)) +  geom_point(aes(shape =Method), size = 2) + geom_hline(yintercept = 0.95)  +   scale_y_continuous(breaks = sort(c(seq( 0.6, 1, length.out=5), 0.95)),limits = c(0.6,1)) + scale_shape_manual(values =c(2,3,4,5,6,7,8)) +   scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash','solid')) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
          labs(title=paste0(figureTitle," (D)") ,
         y = "Coverage",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_npg() + scale_fill_npg()

    
    
    # Figures
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,1],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    
    df = df[,c("CARE_Efron_BIC","MRmix","Egger","RAPS","ConMix","MRLasso","theta")]
    colnames(df) = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso","theta")
    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso"))
    
    
    #df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","dIVW","Mode","RAPS","ConMix","PRESSO"))
    # Power # Bias
    
    #
    p5 =  ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) + geom_line(aes(linetype = Method)) +
    geom_point(aes(shape =Method), size = 2)+geom_hline(yintercept = 0.05) +   scale_y_continuous(breaks = sort(c(seq(0,1, length.out=5), 0.05))) + scale_shape_manual(values =c(2,8,9,10,17,18)) + scale_linetype_manual(values = c('solid', '1F', 'F1', '4C88C488', '12345678','twodash'))  + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
    labs(title=paste0(figureTitle," (E)") ,
         y = "Power",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + guides(col =
         guide_legend(nrow = 1),shape =guide_legend(nrow = 1) )  + scale_color_lancet() + scale_fill_lancet()
    
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(abs(final.out[,"Bias"]),tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    
    df = df[,c("CARE_Efron_BIC","MRmix","Egger","RAPS","ConMix","MRLasso","theta")]
    colnames(df) = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso","theta")
    df[,1:6] = abs(df[,1:6]) #/ df[,7]
    #df = df[df[,"theta"]!=0,]
    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso"))
    
    #colnames(df) = c("ms_MA_AIC","ms_AIC","ms_MA_BIC","rp_MA_AIC","rp_AIC","rp_MA_BIC","theta")

    df.plot = df %>% gather(key = Method, value = Power, -theta)
    df.plot[df.plot[,3]>1,3] = 1

    #df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","dIVW","Mode","RAPS","ConMix","PRESSO"))
    
    # Power # Bias
    
    #stat_smooth(aes(fill=Method, linetype=Method), size=0.5, se = FALSE, level = 0.95) +
    p6 = ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) + geom_line(aes(linetype = Method)) +  geom_point(aes(shape =Method), size = 2) + scale_shape_manual(values =c(2,8,9,10,17,18)) + scale_linetype_manual(values = c('solid', '1F', 'F1', '4C88C488', '12345678','twodash')) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
          labs(title=paste0(figureTitle," (F)") ,
         y = "Absolute Bias",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_lancet() + scale_fill_lancet()

    
    
    
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,"MSE"],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    df = df[,c("CARE_Efron_BIC","MRmix","Egger","RAPS","ConMix","MRLasso","theta")]
    colnames(df) = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso","theta")
    
    #if(min(df[,"MR-Egger"]) > max(df[,c("CARE","MR-mix","RAPS","ContMix","MR-Lasso")]) ) {
    #    df[,"MR-Egger"] = max(df[,c("CARE","MR-mix","RAPS","ContMix","MR-Lasso")]) * 1.2
    #}

    df.plot = df %>% gather(key = Method, value = Power, -theta)
    df.plot[df.plot[,"Power"]>MSE_trunc,"Power"] = MSE_trunc
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso"))
    
    #
    p7 = ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) +geom_line(aes(linetype = Method)) +  geom_point(aes(shape =Method), size = 1.5, alpha=0.85)  + scale_shape_manual(values =c(2,8,9,10,17,18)) + scale_linetype_manual(values = c('solid', '1F', 'F1', '4C88C488', '12345678','twodash')) +  ylim(0,MSE_trunc) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
    axis.text = element_text( size = rel(1.2)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
    legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
    labs(title=paste0(figureTitle," (G)") ,
    y = "Mean Squared Error",
    x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_lancet() + scale_fill_lancet()

    
    
    
    out = matrix(NA, length(files),23)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*", "", tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,"Coverage"],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    df = df[,c("CARE_Efron_BIC","MRmix","Egger","RAPS","ConMix","MRLasso","theta")]
    colnames(df) = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso","theta")

    
    df.plot = df %>% gather(key = Method, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","MR-mix","MR-Egger","RAPS","ContMix","MR-Lasso"))
    
    df.plot[df.plot[,3]<=0.6,3] = 0.6
    
    #scale_y_continuous(breaks = sort(c(seq( 0.6, 1, length.out=5), 0.95))) +
    p8 =  ggplot(df.plot, aes(x=theta, y=Power, color=Method,group = Method)) +geom_line(aes(linetype = Method)) +  geom_point(aes(shape =Method), size = 2) + geom_hline(yintercept = 0.95) +    scale_y_continuous(breaks = sort(c(seq( 0.6, 1, length.out=5), 0.95)),limits = c(0.6,1)) + scale_shape_manual(values =c(2,8,9,10,17,18)) + scale_linetype_manual(values = c('solid', '1F', 'F1', '4C88C488', '12345678','twodash')) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
          labs(title=paste0(figureTitle," (H)") ,
         y = "Coverage",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_lancet() + scale_fill_lancet()
  
    
    g_legend<-function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)}
    mylegend<-g_legend(p1)
    mylegend2<-g_legend(p5)
    
    library(gridExtra)
    
    
    #setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MR-biascorrection/")
    
    jpeg(paste0(saveTitle,".jpeg"),  units="in", width=14, height=8, res=300)
    grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    p4 + theme(legend.position="none"),
    nrow=1),  mylegend, arrangeGrob(p5 + theme(legend.position="none"),
    p6 + theme(legend.position="none"),
    p7 + theme(legend.position="none"),
    p8 + theme(legend.position="none"),
    nrow=1),  mylegend2, nrow=4,heights=c(10, 1.5,10,1.5))
    
    dev.off()
    
}



simSetting = "simRes13"
prop = 0.5
N = "5e+05"
for (M in c(1e+05)) {
    M1 = "1e\\+05"
    files = list.files(getwd())
    files = files[grepl(simSetting,files)]
    files = files[grepl(".RData",files)]
    files = files[grepl(paste0("M", M1),files)]
    #files = files[grepl(paste0("N", N,  "prop_invalid",prop),files)]
    figureTitle = ""
    draw.figure(figureTitle,paste0(simSetting,"_N",N,"_", prop, "_M", M))  
}
for (M in c(1e+05)) {
    M1 = "1e\\+05"
    files = list.files(getwd())
    files = files[grepl(simSetting,files)]
    files = files[grepl(".RData",files)]
    files = files[grepl(paste0("M", M1),files)]
    #files = files[grepl(paste0("N", N,  "prop_invalid",prop),files)]
    figureTitle = ""
    draw.figure(figureTitle,paste0(simSetting,"_N",N,"_", prop, "_M", M))  
}

for (M in c(50000,10000,5000,1000)) {
    files = list.files(getwd())
    files = files[grepl(simSetting,files)]
    files = files[grepl(".RData",files)]
    files = files[grepl(paste0("M", M),files)]
    #files = files[grepl(paste0("N", N,  "prop_invalid",prop),files)]
    figureTitle = ""
    draw.figure(figureTitle,paste0(simSetting,"_N",N,"_", prop, "_M", M))  
}