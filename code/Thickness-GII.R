#!/usr/bin/env Rscript
#Created by N. Crossley, U.C. Chile, modified by A. Zugman, NIMH/NIH
#Meta-analysis and figure creation for 
#"Country-level gender inequality is associated with structural differences in the 
#brains of women and men"
#
#load packages
library(metafor)
library(ggplot2)
library(egg)
library(cowplot)
library(ggseg)
library(here)
library(rstudioapi)

#------------- INPUT NEEDED---------------------
# Change the directory as needed
if ( rstudioapi::isAvailable() ){
  script_path <- dirname( rstudioapi::getActiveDocumentContext()$path )
  rootdir <- dirname(script_path)
  setwd (rootdir)
} else {
  script_path <- paste0(here(),"/","code","/")
  rootdir <- dirname(dirname(script_path))
  setwd (rootdir)
}

load(paste0(rootdir,"/","tables","/","dataformeta.Rdata"))
#for the output
dir.create(paste0(rootdir,"/","figures"))


#----------- Meta-analysis ---------------------
Am=matrix(NA,numTH,2)
for (z in 1:numTH){
  M<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GI)),method = "PM")
  Am[z,1]=A$beta[2]
  Am[z,2]=A$pval[2]
  rm(M,A)
}
rm(z)

# For regional
Am_ROI=Am
print("pval for Hemispheres:")
print(Am_ROI[c(35,70),])
Am_ROI=Am_ROI[-c(35,70),] # excluding the global measures for the FDR correction
Am_fdr=p.adjust(Am_ROI[,2],"fdr")
print("p<0.05: ")
print (Am_fdr[Am_fdr < 0.05])

# Controlling for GDP
Amg=matrix(NA,numTH,2)
for (z in 1:numTH){
  Mg<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  Ag<-rma.uni(yi=Mg$yi,vi=Mg$vi,mods=cbind(as.numeric(GI),log(as.numeric(GDP))),method = "PM")
  Amg[z,1]=Ag$beta[2]
  Amg[z,2]=Ag$pval[2]
  rm(Mg,Ag)
}
rm(z)
Am_ROIg<-Amg
Am_ROIg<-Am_ROIg[-c(35,70),]
Am_fdrg<-p.adjust(Am_ROIg[,2],"fdr")
print("p<0.05 controlling for GDP: ")
print (Am_fdrg[Am_fdrg < 0.05])

#----------------------------------- Figures -----------------------------------
z=35
M1<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A1<-rma.uni(yi=M1$yi,vi=M1$vi,mods=as.numeric(GI),method = "PM")

Tright1<-as.data.frame(cbind(as.numeric(GI),M1$yi,as.numeric(DEM$N)))
Tright1$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 or 400? to improve visibility
Pex_r1<-predict.rma(A1, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp<- ggplot(Tright1, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=15, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Difference in cortical thickness (residuals)") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=1, y=-0.08, label= "paste(italic(P), \"-value = 0.006\")", parse = TRUE,size=7,hjust=1) + 
  annotate("text",x=1, y=-0.095, label= "paste(italic(I) ^ 2, \" = 24.38%\")", parse = TRUE,size=7,hjust=1) + 
  annotate("text",x=1, y=-0.11, label= "paste(italic(R) ^ 2, \" = 3.97%\")", parse = TRUE,size=7,hjust=1) + 
  annotate("label",x=-1.3, y=-0.1, label= "FEMALE > MALE", parse = TRUE,size=7) + 
  annotate("label",x=-1.3, y=0.1, label= "FEMALE < MALE", parse = TRUE,size=7) + 
  ggtitle("A. Right Hemisphere Cortical Thickness") +
  theme(plot.title=element_text(size=22, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=18),  # X axis title
        axis.title.y=element_text(size=18),  # Y axis title
        axis.text.x=element_text(size=15),  # X axis text
        axis.text.y=element_text(size=15),
        legend.title = element_text(size=19), 
        legend.text = element_text(size=16),
        legend.key=element_rect(fill='white'),
        legend.position="bottom",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_r1$pred[1], xend = 1, yend = Pex_r1$pred[2]),size=1,linetype="solid",colour="black") 

ggsave('fig1.pdf',plot=bp,path=paste0(rootdir,"/","figures","/"),dpi=800,width=14,height=7, units="in" )
command = paste0("open ",rootdir,"/","figures","/","fig1.pdf") 
system(command)

#----------------------- ROI figures -------------------------------------------

z=2
M2<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A2<-rma.uni(yi=M2$yi,vi=M2$vi,mods=as.numeric(GI),method = "PM")

Tright2<-as.data.frame(cbind(as.numeric(GI),M2$yi,as.numeric(DEM$N)))
Tright2$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r2<-predict.rma(A2, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_r1<- ggplot(Tright2, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=10, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab(" ") + ylab("Difference in cortical thickness") + theme_light() +  
  annotate("text",x=0.3, y=-0.18, label= "paste(italic(P), \"-value (fdr) = 0.041\")", parse = TRUE,size=5) + 
  annotate("text",x=0.2, y=-0.21, label= "paste(italic(I) ^ 2, \" = 16.89%, \")", parse = TRUE,size=5,hjust=1) + 
  annotate("text",x=0.2, y=-0.21, label= "paste(italic(R) ^ 2, \" = 30.7%\")", parse = TRUE,size=5,hjust=0) + 
  ylim(c(-0.21,0.22))+
  ggtitle("A. RIGHT CAUDAL\n    ANTERIOR CINGULATE") +
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=13),  # X axis title
        axis.title.y=element_text(size=13),  # Y axis title
        axis.text.x=element_text(size=9),  # X axis text
        axis.text.y=element_text(size=9),
        legend.title = element_text(size=13), 
        legend.text = element_text(size=12),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_r2$pred[1], xend = 1, yend = Pex_r2$pred[2]),size=0.5,linetype="solid",colour="black") 
rm(M2,A2)

# right medial orbitofrontal
z=13
M2<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A2<-rma.uni(yi=M2$yi,vi=M2$vi,mods=as.numeric(GI),method = "PM")

Tright3<-as.data.frame(cbind(as.numeric(GI),M2$yi,as.numeric(DEM$N)))
Tright3$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r3<-predict.rma(A2, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_r2<- ggplot(Tright3, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=10, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Difference in cortical thickness") + theme_light() + 
  annotate("text",x=0.3, y=-0.18, label= "paste(italic(P), \"-value (fdr) = 0.013\")", parse = TRUE,size=5) + 
  annotate("text",x=0.3, y=-0.21, label= "paste(italic(I) ^ 2, \" = 20.92%, \")", parse = TRUE,size=5,hjust=1) + 
  annotate("text",x=0.3, y=-0.21, label= "paste(italic(R) ^ 2, \" = 28.43%\")", parse = TRUE,size=5,hjust=0) + 
  ylim(c(-0.21,0.22))+
  ggtitle("B. RIGHT MEDIAL\n    ORBITOFRONTAL") +
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=13),  # X axis title
        axis.title.y=element_text(size=13),  # Y axis title
        axis.text.x=element_text(size=9),  # X axis text
        axis.text.y=element_text(size=9),
        legend.title = element_text(size=13), 
        legend.text = element_text(size=12),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_r3$pred[1], xend = 1, yend = Pex_r3$pred[2]),size=0.5,linetype="solid",colour="black") +
  rm(M2,A2)

# left lateral occipital
z=45
M2<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A2<-rma.uni(yi=M2$yi,vi=M2$vi,mods=as.numeric(GI),method = "PM")

Tright4<-as.data.frame(cbind(as.numeric(GI),M2$yi,as.numeric(DEM$N)))
Tright4$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r4<-predict.rma(A2, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_r3<- ggplot(Tright4, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=10, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab(" ") + ylab("Difference in cortical thickness") + theme_light() + 
  annotate("text",x=0.3, y=-0.18, label= "paste(italic(P), \"-value (fdr) = 0.043\")", parse = TRUE,size=5) + 
  annotate("text",x=0.3, y=-0.21, label= "paste(italic(I) ^ 2, \" = 18.07%, \")", parse = TRUE,size=5,hjust=1) + 
  annotate("text",x=0.3, y=-0.21, label= "paste(italic(R) ^ 2, \" = 31.78%\")", parse = TRUE,size=5,hjust=0) + 
  ylim(c(-0.21,0.22))+
  ggtitle("C. LEFT LATERAL\n    OCCIPITAL") +
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=13),  # X axis title
        axis.title.y=element_text(size=13),  # Y axis title
        axis.text.x=element_text(size=9),  # X axis text
        axis.text.y=element_text(size=9),
        legend.title = element_text(size=13), 
        legend.text = element_text(size=12),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_r4$pred[1], xend = 1, yend = Pex_r4$pred[2]),size=0.5,linetype="solid",colour="black") 
rm(M2,A2)

someData = data.frame(region = c("caudal anterior cingulate","cuneus","fusiform","medial orbitofrontal","paracentral","posterior cingulate","isthmus cingulate","rostral anterior cingulate",
                                 "superior frontal","precuneus","parahippocampal","pericalcarine","temporal pole",
                                 "lingual","superior parietal","isthmus cingulate","superior frontal",
                                 "corpus callosum","entorhinal"),p = c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2), stringsAsFactors = FALSE)
AC=ggseg(.data=someData, atlas=dk, colour="black", size=0.3,mapping=aes(fill=p),hemisphere="right",view="medial")+theme_void()+theme(legend.position="none")+scale_fill_continuous(low="deepskyblue1", high="grey95")

someData2 = data.frame(region = c("medial orbitofrontal","caudal anterior cingulate","cuneus","fusiform","paracentral","posterior cingulate","isthmus cingulate","rostral anterior cingulate",
                                  "superior frontal","precuneus","parahippocampal","pericalcarine",
                                  "lingual","superior parietal","lateral orbitofrontal","temporal pole","isthmus cingulate",
                                  "corpus callosum","entorhinal"),p = c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2), stringsAsFactors = FALSE)
#someData2 = data.frame(region = c("medial orbitofrontal"),  p = 10,  stringsAsFactors = FALSE)
OF=ggseg(.data=someData2, atlas=dk, colour="black", size=0.3,mapping=aes(fill=p),hemisphere="right",view="medial")+theme_void()+theme(legend.position="none") +scale_fill_continuous(low="deepskyblue1", high="grey95")

someData3 = data.frame(region = c("bankssts","caudal middle frontal", 
                                  "lateral occipital","lateral orbitofrontal", "middle temporal",
                                  "fusiform","inferior parietal","inferior temporal",
                                  "pars opercularis","pars orbitalis" ,"pars triangularis",
                                  "postcentral","precentral" ,"rostral middle frontal",
                                  "superior frontal","superior parietal","superior temporal",
                                  "supramarginal","temporal pole","transverse temporal","insula"
),
p=c(2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),stringsAsFactors = FALSE)
LO=ggseg(.data=someData3, atlas=dk, colour="black", size=0.3,mapping=aes(fill=p),hemisphere="left",view="lateral")+theme_void()+theme(legend.position="none") +scale_fill_continuous(low="deepskyblue1", high="grey95")


BR=plot_grid(AC,OF,LO,nrow=1)

title <- ggdraw() + 
  draw_label(
    "Regional differences in cortical thickness",
    fontface = 'bold',x = 0,hjust = 0,size = 18) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
legend <- get_legend(
  # create some space to the left of the legend
  bp_r2 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

A=ggdraw(bp_r1)  +    draw_plot(AC,x=-0.18,y=.28,scale=0.38)
B=ggdraw(bp_r2 + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_text(size=15))) + draw_plot(OF,x=-0.26,y=.28,scale=0.38) 
C=ggdraw(bp_r3 +theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank() )) +draw_plot(LO,x=-0.26,y=.28,scale=0.38)

F3c=plot_grid(A,B,C , align="vh",nrow=1)

F3ct=plot_grid(title,legend,F3c,ncol=1,rel_heights=c(0.1,0.1,1))
ggsave('fig2.pdf',plot=F3ct,path=paste0(rootdir,"/","figures","/"),dpi=800,width=14,height=7, units="in")
command = paste0("open ",rootdir,"/","figures","/","fig2.pdf") 
system(command)
#print (F3ct)
