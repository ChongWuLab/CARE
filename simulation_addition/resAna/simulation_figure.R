library(tidyr)
library(ggplot2)
library(ggsci)

setwd("/rsrch5/home/biostatistics/chongwulab/wzhang24/CARE/simulation_addition/summaryRes/")


draw.figure <- function(figureTitle,saveTitle,MSE_trunc =0.003) {
    #parameter:
    #MSE_trunc = 0.003
    out = matrix(NA, length(files),25)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        #file = "CD2_simRes11theta-0.2prop_invalid0.3.RData"
        load(file)
        final.out
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*","",tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,1],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    
    df = data.frame(out)
    
    df = df[,c("CARE_Efron_BIC","CARE2_Efron_BIC","CARE3_Efron_BIC", "theta")]
    
    colnames(df) = c("CARE","CARE_l1","CARE_l1_2", "theta")
    
    df.plot = df %>% gather(key = eta, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","CARE_l1","CARE_l1_2"))
    # Power # Bias
    
    #
    p1 = ggplot(df.plot, aes(x=theta, y=Power, color=eta, group = eta)) + geom_line(aes(linetype = eta)) +
    geom_point(aes(shape =eta), size = 2)+geom_hline(yintercept = 0.05) +   scale_y_continuous(breaks = sort(c(seq(0,1, length.out=5), 0.05))) + scale_shape_manual(values =c(2,3,4,5,6,7)) +  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash')) + theme_bw() +
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
         
    
    
    out = matrix(NA, length(files),25)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*","",tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(abs(final.out[,"Bias"]),tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    df = data.frame(out)
    df = df[,c("CARE_Efron_BIC","CARE2_Efron_BIC","CARE3_Efron_BIC", "theta")]
    colnames(df) = c("CARE","CARE_l1","CARE_l1_2", "theta")
    df.plot = df %>% gather(key = eta, value = Power, -theta)
    df.plot[df.plot[,3]>1,3] = 1

    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","CARE_l1","CARE_l1_2"))
    # Power # Bias
    
    
    #stat_smooth(aes(fill=Method, linetype=Method), size=0.5, se = FALSE, level = 0.95) +
    p2 = ggplot(df.plot, aes(x=theta, y=Power, color=eta,group = eta)) + geom_line(aes(linetype = eta)) +  geom_point(aes(shape =eta), size = 2) + scale_shape_manual(values =c(2,3,4,5,6,7)) +  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash')) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
          labs(title=paste0(figureTitle," (B)") ,
         y = "Absolute Bias",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_npg() + scale_fill_npg()

    
    out = matrix(NA, length(files),25)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*","",tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,"MSE"],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    df = data.frame(out)
    
    df = df[,c("CARE_Efron_BIC","CARE2_Efron_BIC","CARE3_Efron_BIC", "theta")]
    colnames(df) = c("CARE","CARE_l1","CARE_l1_2", "theta")
  
    df.plot = df %>% gather(key = eta, value = Power, -theta)
    
    df.plot[df.plot[,"Power"]>MSE_trunc,"Power"] = MSE_trunc

    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","CARE_l1","CARE_l1_2"))
    
    #
    p3 = ggplot(df.plot, aes(x=theta, y=Power, color=eta,group = eta)) +geom_line(aes(linetype = eta)) +  geom_point(aes(shape =eta), size = 1.5, alpha=0.85)  + scale_shape_manual(values =c(2,3,4,5,6,7)) +   scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash')) + ylim(0,MSE_trunc) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
    axis.text = element_text( size = rel(1.2)),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
    legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
    labs(title=paste0(figureTitle," (C)") ,
    y = "Mean Squared Error",
    x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_npg() + scale_fill_npg()
    
    
    
    out = matrix(NA, length(files),25)
    
    for(j in 1:length(files)) {
        
        file = files[j]
        load(file)
        
        tmp = gsub("prop.*","",file)
        tmp = gsub(".*theta","",tmp)
        tmp = gsub("N.*","",tmp)
        tmp = as.numeric(tmp)
        out[j,] = c(final.out[,"Coverage"],tmp)
    }
    
    colnames(out) = c(rownames(final.out),"theta")
    df = data.frame(out)
    
    df = df[,c("CARE_Efron_BIC","CARE2_Efron_BIC","CARE3_Efron_BIC", "theta")]
    colnames(df) = c("CARE","CARE_l1","CARE_l1_2", "theta")
  
    df.plot = df %>% gather(key = eta, value = Power, -theta)
    
    df.plot[,2] <- factor( df.plot[,2],levels = c("CARE","CARE_l1","CARE_l1_2"))
    
    df.plot[df.plot[,3]<=0.6,3] = 0.6

    #
    p4 = ggplot(df.plot, aes(x=theta, y=Power, color=eta,group = eta)) +geom_line(aes(linetype = eta)) +  geom_point(aes(shape =eta), size = 2) + geom_hline(yintercept = 0.95)  +   scale_y_continuous(breaks = sort(c(seq( 0.6, 1, length.out=5), 0.95)),limits = c(0.6,1)) + scale_shape_manual(values =c(2,3,4,5,6,7)) +   scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash', 'longdash','twodash')) + theme_bw() +
    theme(axis.title = element_text(face="bold", size = rel(1.2)),
          axis.text = element_text( size = rel(1.2)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0, unit='pt')),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", margin = margin(r = 10, unit = 'pt'), size = rel(1.2)),
          legend.text = element_text(margin = margin(r = 10, unit = 'pt'), size = rel(1.2) ) )  +
          labs(title=paste0(figureTitle," (D)") ,
         y = "Coverage",
         x = expression(paste("Causal effect size (", beta,")",sep=""))) + scale_color_npg() + scale_fill_npg()

    
    
    g_legend<-function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)}
    mylegend<-g_legend(p1)
    
    library(gridExtra)
    
    
    #setwd("/Users/uniquechong/Dropbox (Personal)/FSU_research/Undergoing/MR-biascorrection/")
    
    jpeg(paste0(saveTitle,".jpeg"),  units="in", width=14, height=4, res=300)
    grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    p3 + theme(legend.position="none"),
    p4 + theme(legend.position="none"),
    nrow=1),  mylegend, nrow=2,heights=c(10, 1.5))
    
    dev.off()
    
}




simSetting = "simRes20"
prop = 0.5
    
files = list.files(getwd())

files = files[grepl(simSetting,files)]
files = files[grepl(paste0("prop_invalid",prop),files)]

figureTitle = ""

draw.figure(figureTitle,paste0("Sup_",simSetting,"_",prop),0.001)


