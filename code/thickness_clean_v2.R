library(metafor)
# libraries for figures
library(ggplot2)
library(egg)
library(cowplot)
library(ggseg)

#------------- INPUT NEEDED---------------------
# Change the directory as needed
CP<-' /Users/zugmana2/Documents/meta-gender/FOR_ZUGMAN/'
#-----------------------------------------------

# List of included studies
COHT<-read.csv(paste0(CP,'TABLE_all_used_v2.csv'))

# ------------------------------------------------------------------------------
## The Gender inequality index combined
# Average of the z-scores (flip signal for World Bank)

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
  #print(temp)
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

# Relationship between GDP and GI in the sample
GDPandGI=cor.test(GDP,GI)

#---------------------LOADING DATA---------------------------
## DATA structure
TEST<-read.csv(paste0(CP,COHT$Site[1],'/DATA.csv'))
NTEMP=names(TEST)
TH=4:73 # Thickness specific indices
rm(TEST)

## THICKNESS
# Variables Thickness (means and SD per gender)
mTHm=matrix(NA,length(COHT$Site),length(TH))
sdTHm=matrix(NA,length(COHT$Site),length(TH))
mTHf=matrix(NA,length(COHT$Site),length(TH))
sdTHf=matrix(NA,length(COHT$Site),length(TH))
# Unresidualized (used for the analyses exploring changes in men and separately in women)
mTHmur=matrix(NA,length(COHT$Site),length(TH)) # unresidualized
sdTHmur=matrix(NA,length(COHT$Site),length(TH)) # idem
mTHfur=matrix(NA,length(COHT$Site),length(TH)) # idem
sdTHfur=matrix(NA,length(COHT$Site),length(TH)) # idem

DEMt=matrix(NA,length(COHT$Site),9)
DEMt[,1]=COHT$Site
for(z in 1:length(COHT$Site)){
  DT<-read.csv(paste0(CP,COHT$Site[z],'/DATA.csv'))
  
  ## Cleaning the data
  # Loading matrix with QC measures
  OUTL<-read.csv(paste0(CP,COHT$Site[z],'/DATA_outliers_cortical.csv'))
  
  # Excluding subjects after visual QC
  te=which(DT$QCi>=3)
  if(length(te)>0){
    DT<-DT[-c(te),]
    OUTL<-OUTL[-c(te),]
  }

  # Excluding subjects after automatic QC
  sp<-which(rowSums(OUTL)>0)
  if(length(sp)>0){
    DT<-DT[-sp,]
  }
  rm(OUTL,te,sp)
  
  # alas, the INPD flipped data
  if(is.element(COHT$Site[z],c("INPD_SP","INPD_POA"))){
    im=which(DT$Sex=="f")
    ife=which(DT$Sex=="m")
  }else{
    im=which(DT$Sex=="m")
    ife=which(DT$Sex=="f")}
  
  #---------------Demographics of analysed data-----------
  DEMt[z,2]=length(DT$S)
  DEMt[z,3]=length(im)
  DEMt[z,4]=length(ife)
  # Ages disregard the months 
  DEMt[z,5]=mean(floor(DT$Age))
  DEMt[z,6]=max(floor(DT$Age))
  DEMt[z,7]=min(floor(DT$Age))
  DEMt[z,8]=mean(floor(DT$Age[ife]))
  DEMt[z,9]=mean(floor(DT$Age[im]))
  
  # Regressing out the age
  for (x in 1:length(TH)) {
    a=DT[,TH[x]]
    ch1=which(a==0) # just in case a region failed to be identified
    if(length(ch1)>0){
      print("DANGER")
      TEST=rbind(TEST,c(z,x))
    }
    a[ch1]=NA
    rm(ch1)
    c=floor(DT$Age)
    wisn=setdiff(1:length(a),which(is.na(a))) # otherwise it scrambles the vectors
    lres=lm(formula= a ~ c )
    
    if( ( length(a)-length(wisn) ) >0 ){
      lresiduals=matrix(NA,length(a),1)
      lresiduals[wisn]=lres$residuals
    }else{
      lresiduals=lres$residuals
    }
    
    mTHm[z,x]=mean(lresiduals[im],na.rm=TRUE)
    sdTHm[z,x]=sd(lresiduals[im],na.rm=TRUE)
    mTHf[z,x]=mean(lresiduals[ife],na.rm=TRUE)
    sdTHf[z,x]=sd(lresiduals[ife],na.rm=TRUE)
    
    # Age not regressed out
    mTHmur[z,x]=mean(a[im],na.rm=TRUE)
    sdTHmur[z,x]=sd(a[im],na.rm=TRUE)
    mTHfur[z,x]=mean(a[ife],na.rm=TRUE)
    sdTHfur[z,x]=sd(a[ife],na.rm=TRUE)
    
    rm(a,c,lres,lresiduals,wisn)
  }
  rm(DT,ife,im,x)
}
rm(z)

# Just cleaning up the Demographics and extended Demographics
DEM=as.data.frame(DEMt)
names(DEM)=c("Sample","N","NMale","NFem","AGEm","MaxAge","MinAge","FEM_AGEm","MALE_AGEm")
rm(DEMt)

# Excluding any subject with only males or females after clean-up
toex=union(which(as.numeric(DEM$NFem)<2),which(as.numeric(DEM$NMale)<2))
if(length(toex)>0){
  mTHm<-mTHm[-c(toex),]
  sdTHm<-sdTHm[-c(toex),]
  mTHf<-mTHf[-c(toex),]
  sdTHf<-sdTHf[-c(toex),]
  DEM<-DEM[-c(toex),]
  GI<-GI[-c(toex),]
  GDP<-GDP[-c(toex),]
  COHT<-COHT[-c(toex),]
  mTHmur<-mTHmur[-c(toex),]
  sdTHmur<-sdTHmur[-c(toex),]
  mTHfur<-mTHfur[-c(toex),]
  sdTHfur<-sdTHfur[-c(toex),]
}

#------------------------------------------------------------------------------

# Meta-analysis
Am=matrix(NA,length(TH),2)
for (z in 1:(length(TH))){
  M<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GI)),method = "PM")
  Am[z,1]=A$beta[2]
  Am[z,2]=A$pval[2]
  rm(M,A)
}
rm(z)

# For regional
Am_ROI=Am
Am_ROI=Am_ROI[-c(35,70),] # excluding the global measures for the FDR correction
Am_fdr=p.adjust(Am_ROI[,2],"fdr")

# Controlling for GDP
Amg=matrix(NA,length(TH),2)
for (z in 1:(length(TH))){
  Mg<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  Ag<-rma.uni(yi=Mg$yi,vi=Mg$vi,mods=cbind(as.numeric(GI),log(as.numeric(GDP))),method = "PM")
  Amg[z,1]=Ag$beta[2]
  Amg[z,2]=Ag$pval[2]
  rm(Mg,Ag)
}
rm(z)
Am_ROIg=Amg
Am_ROIg=Am_ROIg[-c(35,70),]
Am_fdrg=p.adjust(Am_ROIg[,2],"fdr")

#===============================================================================
# CHARACTERISING THE SAMPLE
# Demographics
Nr_of_MRI=sum(as.numeric(DEM$N))
Nr_of_femMRI=sum(as.numeric(DEM$NFem))
Nr_of_maleMRI=sum(as.numeric(DEM$NMale))

# Lots of very inefficient code with the figures
COUNTR=unique(COHT$Country)
USAn=which(COHT$Country=="USA")
CHINAn=which(COHT$Country=="China")
EUR=c("Spain","UK","Netherlands","Germany","France","Austria","Italy","Poland","Finland","Sweden","Switzerland")
EURn=vector()
for(x in 1:length(EUR)){
temp1=which(COHT$Country==EUR[x])
EURn=c(EURn,temp1)
}
rm(temp1)
EURn=sort(EURn)
LAT=c("Argentina","Colombia","Cuba","Brazil","Mexico","Chile")
LATn=vector()
for(x in 1:length(LAT)){
  temp1=which(COHT$Country==LAT[x])
  LATn=c(LATn,temp1)
}
rm(temp1)
LATn=sort(LATn)
oHIC=c("Australia","Canada", "Israel","Skorea")
oHICn=vector()
for(x in 1:length(oHIC)){
  temp1=which(COHT$Country==oHIC[x])
  oHICn=c(oHICn,temp1)
}
rm(temp1)
oHICn=sort(oHICn)
oLMIC=c("Russia","India","South Africa","Turkey")
oLMICn=vector()
for(x in 1:length(oLMIC)){
  temp1=which(COHT$Country==oLMIC[x])
  oLMICn=c(oLMICn,temp1)
}
rm(temp1)
oLMICn=sort(oLMICn)

nr_USA=c(sum(as.numeric(DEM$N[USAn])),100*sum(as.numeric(DEM$N[USAn]))/sum(as.numeric(DEM$N)))
nr_CHINA=c(sum(as.numeric(DEM$N[CHINAn])),100*sum(as.numeric(DEM$N[CHINAn]))/sum(as.numeric(DEM$N)))
nr_EUR=c(sum(as.numeric(DEM$N[EURn])),100*sum(as.numeric(DEM$N[EURn]))/sum(as.numeric(DEM$N)))
nr_LAT=c(sum(as.numeric(DEM$N[LATn])),100*sum(as.numeric(DEM$N[LATn]))/sum(as.numeric(DEM$N)))
nr_oHIC=c(sum(as.numeric(DEM$N[oHICn])),100*sum(as.numeric(DEM$N[oHICn]))/sum(as.numeric(DEM$N)))
nr_oLMIC=c(sum(as.numeric(DEM$N[oLMICn])),100*sum(as.numeric(DEM$N[oLMICn]))/sum(as.numeric(DEM$N)))

nCHI=sum(as.numeric(DEM$N[which(COHT$Country=="Chile")]))

TCOUNTR=as.data.frame(rbind(nr_USA,nr_EUR,nr_oHIC,nr_CHINA,nr_LAT,nr_oLMIC),row.names=c("USA","Europe","Other High Income","China","Latin America","Other LMIC"))
names(TCOUNTR)=c("Number of participants","% Total")

nLMIC=sum(TCOUNTR$`Number of participants`[4:6])- nCHI -nr_CHINA[1]
nHIC=sum(TCOUNTR$`Number of participants`[1:3])+nCHI - nr_USA[1]

LHI=as.data.frame(cbind(c(nLMIC,nr_CHINA[1],nHIC,nr_USA[1]),c("","China","","USA"),c("LMIC","LMIC","HIC","HIC")))
names(LHI)=c("N","SAMP","WB")

# Bar chart LMIC/HIC
LHIC=ggplot(LHI,aes(x=WB,y=as.numeric(N),fill=SAMP)) +
  geom_bar(position="stack",stat="identity") + xlab("") +
  ylab("Number of participants") + #coord_flip() +
  ggtitle("C. World Bank Country\nClassification by Income") + labs(fill="Country")+
  scale_fill_viridis_d(labels=c("Other","China","USA"))+
  theme(plot.title=element_text(size=13, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=10),
        legend.title = element_text(size=11), 
        legend.text = element_text(size=11),
        legend.key=element_rect(fill='white'),
        panel.background = element_rect(fill=NA,colour="grey"))+
  scale_fill_discrete(type=c("dimgrey","darkslategray3","darksalmon"),labels=c("Other","China","USA"))

## Bar chart
NLAT=c("Argentina","Brazil","Chile",'Cuba',"Colombia","Mexico")
LAT=which(COHT$Country==NLAT[1])
for(m in 2:length(NLAT)){
      LAT<-union(LAT,which(COHT$Country==NLAT[m])) 
}

NSSA=c("South Africa")
SSA=which(COHT$Country==NSSA[1])

NSA=c("India")
SA=which(COHT$Country==NSA[1])

NME=c("Israel")
ME=which(COHT$Country==NME[1])

NNAM=c("USA","Canada")
NAM=which(COHT$Country==NNAM[1])
for(m in 2:length(NLAT)){
  NAM<-union(NAM,which(COHT$Country==NNAM[m])) 
}

NECA=c("Austria","Spain","Germany","Netherlands","UK","France","Poland","Finland","Sweden","Italy","Switzerland","Russia","Turkey")
ECA=which(COHT$Country==NECA[1])
for(m in 2:length(NECA)){
  ECA<-union(ECA,which(COHT$Country==NECA[m])) 
}

NEAP=c("Australia","China","Skorea")
EAP=which(COHT$Country==NEAP[1])
for(m in 2:length(NEAP)){
  EAP<-union(EAP,which(COHT$Country==NEAP[m])) 
}

B=as.data.frame(cbind(as.numeric(DEM$N),COHT$Country,as.numeric(GI),matrix(NA,length(GI),1)))
names(B)=c("N","Country","GI","Region")
B$Region[EAP]="EA&P"
B$Region[ECA]="Europe and Central Asia"
B$Region[LAT]="Latin America and the Caribb."
B$Region[ME]="ME&NA" 
B$Region[NAM]="NAm" 
B$Region[SA]="SA" 
B$Region[SSA]="SSA"
B$Sample=sample(1:135,replace=FALSE)

B$GIN= (as.numeric(B$GI)+2)*600
  
BG_ALL2=ggplot(B,aes(x=Country)) + geom_bar(aes(fill=Sample,y=as.numeric(N),x=Country),position="stack",stat="identity") + theme(legend.position = "none")+
  geom_point(aes(y=GIN),shape=23,colour="black",fill="white",size=3,stroke=1)+
  facet_grid(~Region,scales="free_x",space="free") +
  scale_colour_grey(start=0.1, end=.9)+
  scale_y_continuous(
    name = "Number of subjects",
    sec.axis = sec_axis(~(./600)-2, name="Gender inequality (z-score)")) +
      theme(axis.text.x = element_text(angle = 90)) +
      ggtitle("A. Studies, countries and gender inequality") + ylab("Number of participants")+
      theme(plot.title=element_text(size=13, 
                                    face="bold", 
                                    hjust=0),  # title
            axis.title.x=element_text(size=12),  # X axis title
            axis.title.y=element_text(size=12),  # Y axis title
            axis.text.x=element_text(size=10),  # X axis text
            axis.text.y=element_text(size=10),
            legend.title = element_text(size=11), 
            legend.text = element_text(size=11),
            legend.key=element_rect(fill='white'),
            panel.background = element_rect(fill=NA,colour="grey")) +
      #scale_fill_grey(start=0.5,end=1)
      #scale_fill_hue(l=60,h=c(1,300),c=70)+
      theme(
        strip.text.x = element_text(
          size = 9, face = "italic"
        )) +
  scale_fill_viridis_c(option="cividis") # magma also looks good

AGEd=as.data.frame(as.numeric(DEM$AGEm))
names(AGEd)="Mean_Age"
MeanAGE=ggplot(AGEd,aes(x=Mean_Age)) + geom_histogram(binwidth=2,color="black", fill="grey80") +
  ggtitle("D. Mean Age and Sex") + xlab("Mean age (years)") + xlim(c(18,40))+
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        legend.title = element_text(size=11), 
        legend.text = element_text(size=10),
        legend.key=element_rect(fill='white'),
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_vline(aes(xintercept=median(AGEd$Mean_Age,na.rm=TRUE)),
             color="black", linetype="dashed", size=0.5)

SEXd=as.data.frame(cbind ( c(sum(as.numeric(DEM$NMale)),sum(as.numeric(DEM$NFem))),c("Men","Women") ))
names(SEXd)=c("N","Sex")
SX=ggplot(SEXd,aes(x=Sex,y=as.numeric(N))) + geom_bar(stat="identity",color="black",fill="gray80")+
  #geom_histogram(binwidth=2,color="black", fill="grey") +
  ggtitle("Sex or\nGender") + ylab("N") + xlab("") +
  theme(plot.title=element_text(size=11, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=8),  # X axis title
        axis.title.y=element_text(size=8),  # Y axis title
        axis.text.x=element_text(size=7),  # X axis text
        axis.text.y=element_text(size=7),
        legend.title = element_text(size=11), 
        legend.text = element_text(size=10),
        legend.key=element_rect(fill='white'),
        panel.background = element_rect(fill=NA,colour="grey")) +
  theme(axis.text.x = element_text(angle = 90))

SX2=ggdraw(MeanAGE)  +    draw_plot(SX,x=0.32,y=-0.02,width=0.9,height=1.3,scale=0.4)

F1=grid.arrange(BG_ALL2,LHIC,SX2, heights=c(1.5,1,1), layout_matrix = rbind(c(1,1,1,1,1),
                                                           #   c(1,1,1,1,1),
                                                       c(NA,NA,NA,NA,2),
                                                       c(NA,NA,NA,NA,3)))

#ggsave('Fig1_new.pdf',plot=F1,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=14,height=12, units="in") 

#===============================================================================
# Reliability
# Checking the reliability of the 3 main results (2,13,35)

rtocheck=c(2,13,35,45) # reliable results
Am_rel=matrix(NA,length(COHT$Site),4)
M1<-escalc("MD",m1i=mTHm[,2],m2i=mTHf[,2],sd1i=sdTHm[,2],sd2i=sdTHf[,2],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
M2<-escalc("MD",m1i=mTHm[,13],m2i=mTHf[,13],sd1i=sdTHm[,13],sd2i=sdTHf[,13],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
M3<-escalc("MD",m1i=mTHm[,45],m2i=mTHf[,45],sd1i=sdTHm[,45],sd2i=sdTHf[,45],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
M4<-escalc("MD",m1i=mTHm[,35],m2i=mTHf[,35],sd1i=sdTHm[,35],sd2i=sdTHf[,35],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
for (z in 1:length(COHT$Site)){
  s=setdiff(1:length(COHT$Site),z)
  A1<-rma.uni(yi=M1$yi[s],vi=M1$vi[s],mods=cbind(as.numeric(GI[s])),method = "PM")
  Am_rel[z,1]=A1$pval[2]
  A2<-rma.uni(yi=M2$yi[s],vi=M2$vi[s],mods=cbind(as.numeric(GI[s])),method = "PM")
  Am_rel[z,2]=A2$pval[2]
  A3<-rma.uni(yi=M3$yi[s],vi=M3$vi[s],mods=cbind(as.numeric(GI[s])),method = "PM")
  Am_rel[z,3]=A3$pval[2]
  A4<-rma.uni(yi=M4$yi[s],vi=M4$vi[s],mods=cbind(as.numeric(GI[s])),method = "PM")
  Am_rel[z,4]=A4$pval[2]
  rm(s,A1,A2,A3,A4)
  print(z)  
}
Am_rel_DF=as.data.frame(Am_rel)
names(Am_rel_DF)=NTEMP[TH[rtocheck]]

# Figures for reliability
REL_35=ggplot(Am_rel_DF,aes(x=rh_MeanThickness_thickness)) + geom_histogram(bins=50,color="black", fill="grey80") +
  ggtitle("A. Right Hemisphere") + xlab("P-value (uncorrected)") + xlim(c(0,0.05))+
  theme(plot.title=element_text(size=12, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=11),  # X axis title
        axis.title.y=element_text(size=11),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_vline(aes(xintercept=Am[35,2]),
             color="black", linetype="dashed", size=0.5)

REL_2=ggplot(Am_rel_DF,aes(x=Am_rel_DF$rh_caudalanteriorcingulate_thickness)) + geom_histogram(bins=50,color="black", fill="grey80") +
  ggtitle("B. Right Caudal\n     Anterior Cingulate") + xlab("P-value (uncorrected)") + xlim(c(0,0.006))+
  theme(plot.title=element_text(size=12, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=11),  # X axis title
        axis.title.y=element_text(size=11),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_vline(aes(xintercept=Am[2,2]),
             color="black", linetype="dashed", size=0.5)

REL_13=ggplot(Am_rel_DF,aes(x=rh_medialorbitofrontal_thickness)) + geom_histogram(bins=50,color="black", fill="grey80") +
  ggtitle("C. Right Medial\n     Orbitofrontal") + xlab("P-value (uncorrected)") + xlim(c(0,0.006))+
  theme(plot.title=element_text(size=12, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=11),  # X axis title
        axis.title.y=element_text(size=11),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_vline(aes(xintercept=Am[13,2]),
             color="black", linetype="dashed", size=0.5)

REL_45=ggplot(Am_rel_DF,aes(lh_lateraloccipital_thickness)) + geom_histogram(bins=50,color="black", fill="grey80") +
  ggtitle("D. Left Lateral\n     Occipital") + xlab("P-value (uncorrected)") + xlim(c(0,0.026))+
  theme(plot.title=element_text(size=12, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=11),  # X axis title
        axis.title.y=element_text(size=11),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_vline(aes(xintercept=Am[45,2]),
             color="black", linetype="dashed", size=0.5)

SF_Rel=grid.arrange(REL_35,REL_2,REL_13,REL_45, layout_matrix = rbind(c(1,1,1,1,2,2,3,3),
                                                                      c(1,1,1,1,NA,4,4,NA)))

#ggsave('SF_Rel_new.pdf',plot=SF_Rel,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=9,height=5, units="in") 


#===============================================================================
# For World Bank indicator only
# World Bank
GIto2=read.csv(paste0(CP,"data_WB.csv"))
t1=which(GIto2$Indicator=="Overall Global Gender Gap Index")
t2=which(GIto2$Subindicator.Type=="Index")
t3=intersect(t1,t2)
GIt2=GIto2[t3,]
rm(t1,t2,t3,GIto2)

GIwb=matrix(NA,length(COHT$Site),1)
for(x in 1:length(GIwb)){
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

  if(temp=="Korea (Republic of)"){
    temp3=which(GIt2$Country.Name=="Korea, Rep.")
  }else {
    temp3=which(GIt2$Country.Name==temp)
  }
  GIwb[x]=GIt2$X2019[temp3]
  rm(temp,temp3)
}

# Meta-analysis
Am_wb=matrix(NA,length(TH),2)
for (z in 1:(length(TH))){
  M<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GIwb)),method = "PM")
  Am_wb[z,1]=A$beta[2]
  Am_wb[z,2]=A$pval[2]
  rm(M,A)
}
Am_ROI_wb=Am_wb
Am_ROI_wb=Am_ROI_wb[-c(35,70),]
Am_fdr_wb=p.adjust(Am_ROI_wb[,2],"fdr")

#===============================================================================
# UNDP
GIt=read.csv(paste0(CP,"2020_statistical_annex_table5_GI.csv"))
GIun=matrix(NA,length(COHT$Site),1)
for(x in 1:length(GIun)){
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
  
  GIun[x]=as.numeric(GIt$GI2019[temp2])
  rm(temp,temp2)
}

# Meta-analysis
Am_un=matrix(NA,length(TH),2)
for (z in 1:(length(TH))){
  M<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  A<-rma.uni(yi=M$yi,vi=M$vi,mods=cbind(as.numeric(GIun)),method = "PM")
  Am_un[z,1]=A$beta[2]
  Am_un[z,2]=A$pval[2]
  rm(M,A)
}
Am_ROI_un=Am_un
Am_ROI_un=Am_ROI_un[-c(35,70),]
Am_fdr_un=p.adjust(Am_ROI_un[,2],"fdr")


#===============================================================================

# Checking that it is a decrease in cortical thickness of females or an increase in males
S=c(2,13,35,45)
DSEX=matrix(NA,length(S),7)
DSEX=as.data.frame(DSEX)
names(DSEX)=c("ROI","Beta Males","P Males","R2 Males","Beta Females", "P Females","R2 Females")
DSEX[,1]=NTEMP[TH[S]]

for(z in 1:length(S)){
  #Males
  mM<-escalc(measure="MN",mi=mTHmur[,S[z]],sdi=sdTHmur[,S[z]],ni=as.numeric(DEM$NMale),vtype="LS")
  mA<-rma.uni(yi=mM$yi,vi=mM$vi,mods=cbind(as.numeric(GI),as.numeric(DEM$MALE_AGEm)),method = "PM")
  
  DSEX[z,2]=mA$beta[2]
  DSEX[z,3]=mA$pval[2]
  DSEX[z,4]=mA$R2

   # Females
   fM<-escalc(measure="MN",mi=mTHfur[,S[z]],sdi=sdTHfur[,S[z]],ni=as.numeric(DEM$NFem),vtype="LS")
   fA<-rma.uni(yi=fM$yi,vi=fM$vi,mods=cbind(as.numeric(GI),as.numeric(DEM$FEM_AGEm)),method = "PM")

   DSEX[z,5]=fA$beta[2]
   DSEX[z,6]=fA$pval[2]
   DSEX[z,7]=fA$R2
  rm(mA,fA,mM,fM)
}


##///////////////////FIGURES//////////////////////////////////

########################################################################
## FIGURE 2 BRIEF REPORT
# Results for global thickness
# Right hemisphere

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
  #geom_segment(aes(x =-1.5, y = Pex_r2$pred[1], xend = 1, yend = Pex_r2$pred[2]),size=0.5,linetype="dotted",colour="black") 

# Male driven or female?
Z=35
mM<-escalc(measure="MN",mi=mTHmur[,Z],sdi=sdTHmur[,Z],ni=as.numeric(DEM$NMale),vtype="LS")
fM<-escalc(measure="MN",mi=mTHfur[,Z],sdi=sdTHfur[,Z],ni=as.numeric(DEM$NFem),vtype="LS")

TMF<-as.data.frame(cbind( c(GI,GI) , c(as.numeric(mM$yi),as.numeric(fM$yi)) , c( as.numeric(DEM$NMale),as.numeric(DEM$NFem)),c(matrix(0,length(DEM$NMale),1), matrix(1,length(DEM$NMale),1))  ))
names(TMF)=c("GI","THICKNESS","N","SEX")
mA<-rma.uni(yi=mM$yi,vi=mM$vi,mods=cbind(as.numeric(GI),as.numeric(DEM$MALE_AGEm)),method = "PM")
fA<-rma.uni(yi=fM$yi,vi=fM$vi,mods=cbind(as.numeric(GI),as.numeric(DEM$FEM_AGEm)),method = "PM")

Pex_m<-predict.rma(mA, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )
Pex_f<-predict.rma(fA, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_mf<- ggplot(TMF, aes(x=GI, y=THICKNESS, size = N,color=factor(SEX))) + 
  geom_point(alpha=0.5) + scale_size_area(max_size=7, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Thickness (mm)") + theme_light() + 
  #dodgerblue4
  scale_colour_manual(values=c("0"="black","1"="dodgerblue4"), name="Sex or gender",labels=c("Male","Female"))+
  #scale_colour_manual(values=c("0"="black","1"="#0066CC"), name="Sex or gender",labels=c("Male","Female"))+
  ggtitle("B. Thickness in women and men") +
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        legend.key=element_rect(fill='white'),
        legend.position="right",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_m$pred[1], xend = 1, yend = Pex_m$pred[2]),size=0.5,linetype="solid",colour="#666666") +
  geom_segment(aes(x =-1.5, y = Pex_f$pred[1], xend = 1, yend = Pex_f$pred[2]),size=0.5,linetype="solid",colour="#0066CC") +
  guides(size="none")


# Variance?

load(file=paste0(CP,"variance_sd_thickness_v3.Rdata"))

z=35
Mv<-escalc("MD",m1i=sdTHm[,z],m2i=sdTHf[,z],sd1i=vsdTHm[,z],sd2i=vsdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
Av<-rma.uni(yi=Mv$yi,vi=Mv$vi,mods=as.numeric(GI),method = "PM")

Trightv<-as.data.frame(cbind(as.numeric(GI),Mv$yi,as.numeric(DEM$N)))
Trightv$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_rv<-predict.rma(Av, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_var<- ggplot(Trightv, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=7, name="Number of\nparticipants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Differences in stand dev of cortical thickness") + theme_light() + theme(legend.position="right") + 
   annotate("label",x=-1.1, y=-0.1, label= "FEMALE > MALE", parse = TRUE) + 
  annotate("label",x=-1.1, y=0.1, label= "FEMALE < MALE", parse = TRUE) + 
   ggtitle("C. Variance difference") +
  theme(plot.title=element_text(size=15, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        legend.key=element_rect(fill='white'),
        legend.position="right",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_rv$pred[1], xend = 1, yend = Pex_rv$pred[2]),size=0.5,linetype="solid",colour="black")
 

F2=grid.arrange(bp,bp_mf,bp_var, widths=c(1,1,1,1,0.3),layout_matrix = rbind(c(1,1,1,2,2),
                                                       c(1,1,1,2,2),
                                                       c(1,1,1,3,3),
                                                       c(1,1,1,3,3)))

#ggsave('Fig2_new.pdf',plot=F2,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=14,height=8, units="in") 



########################################################################
## FIGURE 3 BRIEF REPORT
# Results for ROI thickness
# Right caudal anterior
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

#ggsave('Fig3_new.pdf',plot=F3ct,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=14,height=7, units="in") 




######################################################################################
# SUPPLEMENTARY FIGURE VARIANCE 

z=2
Mv2<-escalc("MD",m1i=sdTHm[,z],m2i=sdTHf[,z],sd1i=vsdTHm[,z],sd2i=vsdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
Av2<-rma.uni(yi=Mv2$yi,vi=Mv2$vi,mods=as.numeric(GI),method = "PM")

Trightv2<-as.data.frame(cbind(as.numeric(GI),Mv2$yi,as.numeric(DEM$N)))
Trightv2$V3[which(as.numeric(DEM$N)>400)]<-400 
Pex_rv2<-predict.rma(Av2, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_var2<- ggplot(Trightv2, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=7, name="Number of\nparticipants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Differences in stand dev thickness") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=0.4, y=-0.2, label= "paste(italic(P), \"-value = 0.65\")", parse = TRUE) + 
  annotate("label",x=-1.3, y=-0.2, label= "FEMALE > MALE", parse = TRUE) + 
  annotate("label",x=-1.3, y=0.2, label= "FEMALE < MALE", parse = TRUE) + 
  ggtitle("Right anterior caudal cingulate") +
  theme(plot.title=element_text(size=13, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_rv2$pred[1], xend = 1, yend = Pex_rv2$pred[2]),size=0.5,linetype="solid",colour="black")

z=13
Mv13<-escalc("MD",m1i=sdTHm[,z],m2i=sdTHf[,z],sd1i=vsdTHm[,z],sd2i=vsdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
Av13<-rma.uni(yi=Mv13$yi,vi=Mv13$vi,mods=as.numeric(GI),method = "PM")

Trightv13<-as.data.frame(cbind(as.numeric(GI),Mv13$yi,as.numeric(DEM$N)))
Trightv13$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_rv13<-predict.rma(Av13, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_var13<- ggplot(Trightv13, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=7, name="Number of\nparticipants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Differences in stand dev thickness") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=0.4, y=-0.09, label= "paste(italic(P), \"-value = 0.91\")", parse = TRUE) + 
  annotate("label",x=-1.25, y=-0.095, label= "FEMALE > MALE", parse = TRUE) + 
  annotate("label",x=-1.25, y=0.13, label= "FEMALE < MALE", parse = TRUE) + 
  ggtitle("Right medial orbitofrontal") +
  theme(plot.title=element_text(size=13, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_rv13$pred[1], xend = 1, yend = Pex_rv13$pred[2]),size=0.5,linetype="solid",colour="black")

z=45
Mv45<-escalc("MD",m1i=sdTHm[,z],m2i=sdTHf[,z],sd1i=vsdTHm[,z],sd2i=vsdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
Av45<-rma.uni(yi=Mv45$yi,vi=Mv45$vi,mods=as.numeric(GI),method = "PM")

Trightv45<-as.data.frame(cbind(as.numeric(GI),Mv45$yi,as.numeric(DEM$N)))
Trightv45$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_rv45<-predict.rma(Av45, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bp_var45<- ggplot(Trightv45, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=7, name="Number of\nparticipants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Differences in stand dev thickness") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=0.3, y=-0.09, label= "paste(italic(P), \"-value = 0.78\")", parse = TRUE) + 
  annotate("label",x=-1.3, y=-0.09, label= "FEMALE > MALE", parse = TRUE) + 
  annotate("label",x=-1.3, y=0.12, label= "FEMALE < MALE", parse = TRUE) + 
  ggtitle("Left lateral occipital") +
  theme(plot.title=element_text(size=13, 
                                face="bold", 
                                hjust=0),  # title
        axis.title.x=element_text(size=12),  # X axis title
        axis.title.y=element_text(size=12),  # Y axis title
        axis.text.x=element_text(size=8),  # X axis text
        axis.text.y=element_text(size=8),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=11),
        legend.key=element_rect(fill='white'),
        legend.position="none",
        panel.background = element_rect(fill=NA,colour="grey")) +
  geom_segment(aes(x =-1.5, y = Pex_rv45$pred[1], xend = 1, yend = Pex_rv45$pred[2]),size=0.5,linetype="solid",colour="black")

title_s1 <- ggdraw() +   draw_label(
  "Difference in standard deviation",
  fontface = 'bold',x = 0,hjust = 0,size = 15) +
  theme(
    plot.margin = margin(0, 0, 0, 10)
  )
legend_s1 <- get_legend(
  # create some space to the left of the legend
  bp_var45 +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

FS_var=plot_grid(title_s1,bp_var2,bp_var13,bp_var45,legend_s1,ncol=1,rel_heights=c(0.1,1,1,1,0.1))

BR2=plot_grid(AC,OF,LO,nrow=3)

FS_all=plot_grid(BR2,FS_var,rel_widths = c(0.3,1),nrow=1)

# ggsave('Fig_SROI_varth_new.pdf',plot=FS_all,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=8,height=10.5, units="in") 

###############################################################3
# Suppl Figure for WB and UNDP
# For right hemisphere
library(cowplot)
z=35
M_un<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A_un<-rma.uni(yi=M_un$yi,vi=M_un$vi,mods=as.numeric(GIun),method = "PM")
A_ung<-rma.uni(yi=M_un$yi,vi=M_un$vi,mods=cbind(as.numeric(GIun),log(as.numeric(GDP))),method = "PM")

Tright_un<-as.data.frame(cbind(as.numeric(GIun),M_un$yi,as.numeric(DEM$N)))
Tright_un$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r_un<-predict.rma(A_un, newmods= rbind(cbind(0.03,25),cbind(0.5,25)) )

bp_un<- ggplot(Tright_un, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=10, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality Index") + ylab("Difference in cortical thickness") + theme_light() +  
  annotate("text",x=0.4, y=-0.09, label= "paste(italic(P), \"-value = 0.026\")", parse = TRUE,size=4) + 
  annotate("text",x=0.4, y=-0.105, label= "paste(italic(I) ^ 2, \" = 22.27%, \")", parse = TRUE,size=4,hjust=1) + 
  annotate("text",x=0.4, y=-0.105, label= "paste(italic(R) ^ 2, \" = 14.27%\")", parse = TRUE,size=4,hjust=0) + 
  #annotate("text",x=0.4, y=-0.1, label= "paste(italic(P), \"-value (GDP corrected) = 0.036\")", parse = TRUE,size=4) + 
  #annotate("text",x=0.4, y=-0.105, label= "paste(italic(I) ^ 2, \" = 22.33%\")", parse = TRUE,size=4) + 
  annotate("label",x=0.08, y=-0.1, label= "FEMALE > MALE", parse = TRUE,size=3) + 
  annotate("label",x=0.08, y=0.1, label= "FEMALE < MALE", parse = TRUE,size=3) + 
  ylim(c(-0.11,0.105))+xlim(c(0,0.5)) +
  ggtitle("UNITED NATIONS") +
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
  geom_segment(aes(x =0.03, y = Pex_r_un$pred[1], xend = 0.5, yend = Pex_r_un$pred[2]),size=0.5,linetype="solid",colour="black")

z=35
M_wb<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A_wb<-rma.uni(yi=M_wb$yi,vi=M_wb$vi,mods=as.numeric(GIwb),method = "PM")
A_wbg<-rma.uni(yi=M_wb$yi,vi=M_wb$vi,mods=cbind(as.numeric(GIwb),log(as.numeric(GDP))),method = "PM")

Tright_wb<-as.data.frame(cbind(as.numeric(GIwb),M_wb$yi,as.numeric(DEM$N)))
Tright_wb$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r_wb<-predict.rma(A_wb, newmods= rbind(cbind(0.64,25),cbind(0.83,25)) )
bp_wb<- ggplot(Tright_wb, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=10, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab("Gender Gap Index") + ylab("Difference in cortical thickness") + theme_light() +  
  annotate("text",x=0.68, y=-0.09, label= "paste(italic(P), \"-value = 0.014\")", parse = TRUE,size=4) + 
  annotate("text",x=0.68, y=-0.105, label= "paste(italic(I) ^ 2, \" = 23.88%, \")", parse = TRUE,size=4,hjust=1) + 
  annotate("text",x=0.68, y=-0.105, label= "paste(italic(R) ^ 2, \" = 6.85%\")", parse = TRUE,size=4,hjust=0) + 
  #annotate("text",x=0.4, y=-0.1, label= "paste(italic(P), \"-value (GDP corrected) = 0.027\")", parse = TRUE,size=4) + 
  #annotate("text",x=0.68, y=-0.105, label= "paste(italic(I) ^ 2, \" = 25.52%\")", parse = TRUE,size=4) + 
  annotate("label",x=0.80, y=-0.095, label= "FEMALE > MALE",parse=TRUE,size=3) + 
  annotate("label",x=0.80, y=0.1, label= "MALE > FEMALE",parse=TRUE,size=3) + 
  ylim(c(-0.11,0.105))+xlim(c(0.61,0.84)) +
  ggtitle("WORLD BANK") + scale_x_reverse()+
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
  geom_segment(aes(x =0.64, y = Pex_r_wb$pred[1], xend = 0.83, yend = Pex_r_wb$pred[2]),size=0.5,linetype="solid",colour="black")


title_s4 <- ggdraw() +   draw_label(
  "Thickness right hemisphere",
  fontface = 'bold',x = 0,hjust = 0,size = 15) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )
legend_s4 <- get_legend(
  # create some space to the left of the legend
  bp_wb +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

FS_both=plot_grid(bp_un,bp_wb,nrow=1)
FS_is=plot_grid(title_s4,FS_both,legend_s4,ncol=1,rel_heights=c(0.1,1,0.1))

# ggsave('Fig_SBOTH_new.pdf',plot=FS_is,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=10,height=5, units="in") 


#////////////////////////////////////////
# Suppl figure results left hemisphere
z=70
M1<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
A1<-rma.uni(yi=M1$yi,vi=M1$vi,mods=as.numeric(GI),method = "PM")

Tright1<-as.data.frame(cbind(as.numeric(GI),M1$yi,as.numeric(DEM$N)))
Tright1$V3[which(as.numeric(DEM$N)>400)]<-400 # capping the N to 200 to improve visibility
Pex_r1<-predict.rma(A1, newmods= rbind(cbind(-1.5,25),cbind(1,25)) )

bpleft<- ggplot(Tright1, aes(x=V1, y=V2, size = V3)) +
  geom_hline(colour="black", yintercept=0,size=0.6,linetype="dotted") +  
  geom_point(alpha=0.3) + scale_size_area(max_size=15, name="Number of participants",breaks=c(100,200,300,400)) +
  xlab("Gender Inequality (z-scores)") + ylab("Difference in cortical thickness") + theme_light() + theme(legend.position="right") + 
  annotate("text",x=0.3, y=-0.09, label= "paste(italic(P), \"-value = 0.10\")", parse = TRUE,size=7) + 
  annotate("text",x=0.3, y=-0.11, label= "paste(italic(I) ^ 2, \" = 22.10%\")", parse = TRUE,size=7) + 
  annotate("label",x=-1.3, y=-0.1, label= "FEMALE > MALE", parse = TRUE,size=7) + 
  annotate("label",x=-1.3, y=0.1, label= "FEMALE < MALE", parse = TRUE,size=7) + 
  ggtitle("Left Hemisphere Cortical Thickness") +
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

#ggsave('Fig_S_left_new.pdf',plot=bpleft,path='/Users/ncrossley/Documents/Laburo/Social_Determinants/Draft_INEQ/Figures/',dpi=800,width=12,height=9, units="in") 

