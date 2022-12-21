library(metafor)
library(ggplot2)
library(egg)
library(cowplot)
library(ggseg)

#------------- INPUT NEEDED---------------------
# Change the directory as needed
CP<-'/Users/ncrossley/Documents/Laburo/Social_Determinants/FREESURFER_2/'
#-----------------------------------------------

# List of studies
COHT<-read.csv(paste0(CP,'TABLE_all_used_v2.csv'))

# UNDP

GIt=read.csv(paste0(CP,"2020_statistical_annex_table5_GI.csv"))
temp1=as.numeric(GIt$GI2019)
zundp=(temp1-mean(temp1,na.rm=TRUE))/sd(temp1,na.rm=TRUE)
rm(temp1)

# World Bank
GIto2=read.csv(paste0(CP,"data_WB.csv"))
t1=which(GIto2$Indicator=="Overall Global Gender Gap Index")
t2=which(GIto2$Subindicator.Type=="Index")
t3=intersect(t1,t2)
GIt2=GIto2[t3,]
temp1=GIt2$X2019
zwb=(temp1-mean(temp1,na.rm=TRUE))/sd(temp1,na.rm=TRUE)
zwb_s=zwb*-1
rm(temp1,t1,t2,t3,GIto2)

GI=matrix(NA,length(COHT$Site))

for(x in 1:length(GI)){
  temp=COHT$Country[x]
  if (temp=="USA"){
    temp="United States"
  } else if(temp=="UK"){
    temp="United Kingdom"
  }else if(temp=="Russia"){
    temp="Russian Federation"
  }else if(temp=="SKorea"){
    temp="Korea (Republic of)"
  }else if(temp=="Skorea"){
    temp="Korea (Republic of)"
  }
  
  temp2=which(GIt$COUNTRY==temp)
  if(temp=="Korea (Republic of)"){
    temp3=which(GIt2$Country.Name=="Korea, Rep.")
  }else {
    temp3=which(GIt2$Country.Name==temp)
  }
  
  GI[x]=(zundp[temp2]+zwb_s[temp3])/2
  rm(temp, temp2, temp3)
}
rm(zundp,zwb,zwb_s,GIt,GIt2)

## GDP per capita
# World Bank
GDPo=read.csv(paste0(CP,"API_NY.GDP.PCAP.CD_DS2_en_csv_v2_3930492.csv"))
GDP=matrix(NA,length(GI),1)
for(x in 1:length(GI)){
  temp=COHT$Country[x]
  if (temp=="USA"){
    temp="United States"
  } else if(temp=="UK"){
    temp="United Kingdom"
  }else if(temp=="Russia"){
    temp="Russian Federation"
  }else if(temp=="SKorea"){
    temp="Korea, Rep."
  }else if(temp=="Skorea"){
    temp="Korea, Rep."
  }
  temp2=which(GDPo$Country.Name==temp)
  GDP[x]=GDPo$X2019[temp2]
  rm(temp,temp2)
}
rm(x,GDPo)

## DATA structure
TEST<-read.csv(paste0(CP,COHT$Site[1],'/DATA.csv'))
NTEMP=names(TEST)
AR=74:143
rm(TEST)

## AREA
# Variables Area
mARm=matrix(NA,length(COHT$Site),length(AR))
sdARm=matrix(NA,length(COHT$Site),length(AR))
mARf=matrix(NA,length(COHT$Site),length(AR))
sdARf=matrix(NA,length(COHT$Site),length(AR))
ARV=matrix(NA,length(COHT$Site),3)

DEMt=matrix(NA,length(COHT$Site),7)
DEMt[,1]=COHT$Site
# Analyses do not include the Colombia and Stellenbosch (2 out of 3 scanners) since they did not send eTIV
for(z in 1:length(COHT$Site)){
  if(is.element(COHT$Site[z],c("Colombia","SouthAfrica_duPlessis_S2","SouthAfrica_duPlessis_S3"))==FALSE){
    DT<-read.csv(paste0(CP,COHT$Site[z],'/DATA.csv'))
    OUTL<-read.csv(paste0(CP,COHT$Site[z],'/DATA_outliers_cortical.csv'))
    te=which(DT$QCi>=3)
    if(length(te)>0){
      DT<-DT[-c(te),]
      OUTL<-OUTL[-c(te),]
    }
    
    sp<-which(rowSums(OUTL)>0)
    if(length(sp)>0){
      DT<-DT[-sp,]
    }
    
    # alas, the INPD flipped data
    if(is.element(COHT$Site[z],c("INPD_SP","INPD_POA"))){
      im=which(DT$Sex=="f")
      ife=which(DT$Sex=="m")
    }else{
      im=which(DT$Sex=="m")
      ife=which(DT$Sex=="f")}
    
    DEMt[z,2]=length(DT$S)
    DEMt[z,3]=length(im)
    DEMt[z,4]=length(ife)
    DEMt[z,5]=mean(floor(DT$Age))
    DEMt[z,6]=max(floor(DT$Age))
    DEMt[z,7]=min(floor(DT$Age))
    
    ARV[z,1]=cor((DT$rh_WhiteSurfArea_area + DT$lh_WhiteSurfArea_area),DT$eTIV)
    ARV[z,2]=mean(DT$eTIV)
    ARV[z,3]=mean(DT$rh_WhiteSurfArea_area + DT$lh_WhiteSurfArea_area) 
    
    for (x in 1:length(AR)) {
      a=DT[,AR[x]]
      ch1=which(a==0)
      if(length(ch1)>0){
        print("DANGER")
        TEST=rbind(TEST,c(z,x))
      }
      a[ch1]=NA
      rm(ch1)
      c=floor(DT$Age)
      c2=floor(DT$Age)^2 # Does not make much difference adding a quadratic term
      d=DT$eTIV
     wisn=setdiff(1:length(a),which(is.na(a)))  
     lres=lm(formula= a ~ c +d ) # regressing out effect of age within the sample, as well as eTIV

      if(length(wisn)>0){
        lresiduals=matrix(NA,length(a),1)
        lresiduals[wisn]=lres$residuals
      }else{
        lresiduals=lres$residuals
      }
      
      mARm[z,x]=mean(lresiduals[im],na.rm=TRUE)
      sdARm[z,x]=sd(lresiduals[im],na.rm=TRUE)
      mARf[z,x]=mean(lresiduals[ife],na.rm=TRUE)
      sdARf[z,x]=sd(lresiduals[ife],na.rm=TRUE)
      rm(pr,nd)
    }
  }
}
DEM=as.data.frame(DEMt)
names(DEM)=c("Sample","N","NMale","NFem","AGEm","MaxAge","MinAge")

# Excluding any subject with only males or females after 
toex=union(which(as.numeric(DEM$NFem)<2),which(as.numeric(DEM$NMale)<2))
mARm<-mARm[-c(toex),]
sdARm<-sdARm[-c(toex),]
mARf<-mARf[-c(toex),]
sdARf<-sdARf[-c(toex),]
DEM<-DEM[-c(toex),]
GI<-GI[-c(toex),]
GDP<-GDP[-c(toex),]
COHT<-COHT[-c(toex),]

# Meta-analysis
Am=matrix(NA,length(AR),2)
for (z in 1:(length(AR))){
  M<-escalc("MD",m1i=mARm[,z],m2i=mARf[,z],sd1i=sdARm[,z],sd2i=sdARf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  A<-rma.uni(yi=M$yi,vi=M$vi,mods=as.numeric(GI),method = "PM", control=list(tau2.max=5000000))
  Am[z,1]=A$beta[2]
  Am[z,2]=A$pval[2]
  rm(M,A)
}
Am_ROI=Am
Am_ROI=Am_ROI[-c(35,70,71),]
Am_fdr=p.adjust(Am_ROI[,2],"fdr")

Amg=matrix(NA,length(AR),2)
for (z in 1:(length(AR))){
  M<-escalc("MD",m1i=mARm[,z],m2i=mARf[,z],sd1i=sdARm[,z],sd2i=sdARf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GI),log(as.numeric(GDP))),method = "PM", control=list(tau2.max=500000))
  Amg[z,1]=A$beta[2]
  Amg[z,2]=A$pval[2]
  rm(M,A)
}
Am_ROIg=Amg
Am_ROIg=Am_ROIg[-c(35,70,71),]
Am_fdrg=p.adjust(Am_ROIg[,2],"fdr")

## Supplementary Figure

# Left hemisphere
z=70
M<-escalc("MD",m1i=mARm[,z],m2i=mARf[,z],sd1i=sdARm[,z],sd2i=sdARf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GI)),method = "PM",control=list(tau2.max=500000))

Tright<-as.data.frame(cbind(as.numeric(GI),M$yi,as.numeric(DEM$N)))
Tright$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r<-predict.rma(A, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp1<- ggplot(Tright, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.7,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=8, name="Number of\nparticipants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Difference in cortical surface area") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=0.3, y=-2100, label= "paste(italic(P), \"-value = 0.38\")", parse = TRUE,size=5) + 
  annotate("text",x=0.3, y=-2700, label= "paste(italic(I) ^ 2, \" = 15.5%\")", parse = TRUE,size=5) + 
  annotate("label",x=-1.3, y=-2500, label= "FEMALE > MALE", parse = TRUE,size=5) + 
  annotate("label",x=-1.3, y=4500, label= "FEMALE < MALE", parse = TRUE,size=5) + 
  ggtitle("Left Hemisphere") +
  theme(plot.title=element_text(size=18, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=15),  # X axis title
        axis.title.y=element_text(size=15),  # Y axis title
        axis.text.x=element_text(size=12),  # X axis text
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_r$pred[1], xend = 1, yend = Pex_r$pred[2]),size=0.5,linetype="solid",colour="black")
rm(M,A)

# Right hemisphere
z=35
M<-escalc("MD",m1i=mARm[,z],m2i=mARf[,z],sd1i=sdARm[,z],sd2i=sdARf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GI)),method = "PM",control=list(tau2.max=500000))

Tright2<-as.data.frame(cbind(as.numeric(GI),M$yi,as.numeric(DEM$N)))
Tright2$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r2<-predict.rma(A, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp2<- ggplot(Tright2, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.7,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=8, name="Number of\nparticipants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Difference in cortical surface area") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=0.3, y=-2100, label= "paste(italic(P), \"-value = 0.40\")", parse = TRUE,size=5) + 
  annotate("text",x=0.3, y=-2700, label= "paste(italic(I) ^ 2, \" = 8.19%\")", parse = TRUE,size=5) + 
  annotate("label",x=-1.3, y=-2500, label= "FEMALE > MALE", parse = TRUE,size=5) + 
  annotate("label",x=-1.3, y=4500, label= "FEMALE < MALE", parse = TRUE,size=5) + 
  ggtitle("Right Hemisphere") +
  theme(plot.title=element_text(size=18, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=15),  # X axis title
        axis.title.y=element_text(size=15),  # Y axis title
        axis.text.x=element_text(size=12),  # X axis text
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_r2$pred[1], xend = 1, yend = Pex_r2$pred[2]),size=0.5,linetype="solid",colour="black")
rm(M,A)

title <- ggdraw() + 
  draw_label(
    "Surface Area Analysis",
    fontface = 'bold',x = 0,hjust = 0,size = 18) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
legend <- get_legend(
  bp2 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "right")
)

A=plot_grid(bp1,bp2,legend,nrow=1,rel_widths=c(1,1,0.3))
B=plot_grid(title,A,ncol=1,rel_heights=c(0.1,1))
# ggsave('Fig_S_surface_new.pdf',plot=B,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures',dpi=800,width=12,height=6, units="in") 

