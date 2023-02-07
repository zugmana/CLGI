#This script will get the data for GI and organize for further use.
#Created by N. Crossley, U.C. Chile, modified by A. Zugman, NIMH/NIH
# Not all data used in the paper can be shared due to individual countries regulation.
##"Country-level gender inequality is associated with structural differences in the 
#brains of women and men"
#
#Library to set wd
library(rstudioapi)
library(metafor)
#------------- INPUT NEEDED---------------------
# Change the directory as needed
script_path <- dirname( rstudioapi::getActiveDocumentContext()$path )
setwd(script_path)
setwd("..")
rootdir <- getwd()
#-----------------------------------------------

# List of included studies
COHT<-read.csv(paste0(rootdir,'/','tables','/','TABLE_all_used_v3.csv'))

# ------------------------------------------------------------------------------
## The Gender inequality index combined
# Average of the z-scores (flip signal for World Bank)

# UNDP
GIt=read.csv(paste0(rootdir,'/','tables','/',"2020_statistical_annex_table5_GI.csv"))
temp1=as.numeric(GIt$GI2019)
zundp=(temp1-mean(temp1,na.rm=TRUE))/sd(temp1,na.rm=TRUE)
rm(temp1)

# World Bank
GIto2=read.csv(paste0(rootdir,'/','tables','/',"data_WB.csv"))
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
GDPo=read.csv(paste0(rootdir,'/','tables','/',"API_NY.GDP.PCAP.CD_DS2_en_csv_v2_3930492.csv"))
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
#
#Now preparing Freesurfer and sex/age
#
#---------------------LOADING DATA---------------------------
## DATA structure
TEST<-read.csv(paste0(rootdir,'/','freesurfertables','/',COHT$Site[1],'/DATA.csv'))
NTEMP=names(TEST)
TH=4:73 # Thickness specific indices
rm(TEST)

## THICKNESS
#Initialize variables
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
  DT<-read.csv(paste0(rootdir,'/','freesurfertables','/',COHT$Site[z],'/DATA.csv'))
  
  ## Cleaning the data
  # Loading matrix with QC measures
  OUTL<-read.csv(paste0(rootdir,'/','freesurfertables','/',COHT$Site[z],'/DATA_outliers_cortical.csv'))
  
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
  
  im=which(DT$Sex=="m")
  ife=which(DT$Sex=="f")
  
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

numTH <- (length(TH))

#Use escalc to get data that will be used in next script
#This gets the effect size for the meta-regression on the next script.
Y <- matrix(NA,length(DEM[,1]), numTH)
Vi <- matrix(NA, length(DEM[,1]), numTH)
for (z in 1:numTH){
  M<-escalc("MD",m1i=mTHm[,z],m2i=mTHf[,z],sd1i=sdTHm[,z],sd2i=sdTHf[,z],n1i=as.numeric(DEM$NMale),n2i=as.numeric(DEM$NFem),vtype="LS")
  Y[ ,z]<-as.numeric(M$yi)
  Vi[ ,z]<-as.numeric(M$vi)
  rm(M)
  }


save(list = c("Y","Vi", "DEM", "GI", "GDP",
              "COHT","numTH"),
     file = paste0(rootdir,"/","tables","/","dataformeta.Rdata"))

