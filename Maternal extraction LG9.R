Ch2=read.table("ChLG9_20000000_21459655.gff3",header = FALSE)

head(Ch2)
Ch2pos=subset(Ch2, V7=="+")

PT=subset(Ch2pos, V3=="processed_transcript")
CE=subset(Ch2pos, V3=="coding_exon")


nreps=length(PT[,1])
start= numeric (nreps)
stop= numeric (nreps)
gene= numeric (nreps)

nreps2=length(CE[,1])
startCE= numeric (nreps2)
stopCE= numeric (nreps2)
exonCE= numeric (nreps2)
geneCE=numeric (nreps2)

for (i in 1:nreps) {
PT1=PT[i,]
start[i]=PT1$V4
stop[i]=PT1$V5	
gene[i]=i
}


k=1
num=0
i=1

for (i in 1:nreps2) {
CE1=subset(CE,V4>=start[i] & V5<=stop[i])
geneCCE=gene[i]
nreps3=length(CE1[,1])

for (k in 1:nreps3) {
num=num+1
startCE[num]=CE1$V4[k]
stopCE[num]=CE1$V5[k]
exonCE[num]=k
geneCE[num]=geneCCE
print(num)
}
}

xxx=cbind(startCE, stopCE, exonCE, geneCE)



Emb6hr=subset(Ch2, V2=="6hr_embryonic" & V3=="region")
EarlyAdultFemale =subset(Ch2, V2=="EarlyAdultFemale" & V3=="region")
EarlyAdultMale =subset(Ch2, V2=="EarlyAdultMale" & V3=="region")
EarlyLastLarval =subset(Ch2, V2=="EarlyLastLarval" & V3=="region")
EarlyPupalFemale =subset(Ch2, V2=="EarlyPupalFemale" & V3=="region")
EarlyPupalMale =subset(Ch2, V2=="EarlyPupalMale"& V3=="region")
LateLastLarval =subset(Ch2, V2=="LateLastLarval"& V3=="region")
MidLastLarval =subset(Ch2, V2=="MidLastLarval"& V3=="region")
MIdPupalFemale =subset(Ch2, V2=="MIdPupalFemale"& V3=="region")
MidPupalMale =subset(Ch2, V2=="MidPupalMale"& V3=="region")



nreps22=length(xxx[,1])

h=1
	
Emb6hr1 =subset(Emb6hr,V4>=startCE[h] & V5<=stopCE[h])
EarlyAdultFemale1=subset(EarlyAdultFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyAdultMale1 =subset(EarlyAdultMale,V4>=startCE[h] & V5<=stopCE[h])
EarlyLastLarval1 =subset(EarlyLastLarval,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalFemale1 =subset(EarlyPupalFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalMale1=subset(EarlyPupalMale,V4>=startCE[h] & V5<=stopCE[h])
LateLastLarval1 =subset(LateLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MidLastLarval1 =subset(MidLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MIdPupalFemale1 =rep(NA,length(Emb6hr1[,1]))
MidPupalMale1 =subset(MidPupalMale,V4>=startCE[h] & V5<=stopCE[h])


what=cbind(as.numeric(paste(Emb6hr1$V6)),as.numeric(paste(EarlyAdultFemale1$V6)),as.numeric(paste(EarlyAdultMale1$V6)),as.numeric(paste(EarlyLastLarval1$V6)),as.numeric(paste(EarlyPupalFemale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(LateLastLarval1$V6)),as.numeric(paste(MidLastLarval1$V6)),as.numeric(paste(MIdPupalFemale1)),as.numeric(paste(MidPupalMale1$V6)), startCE[h], stopCE[h], exonCE[h], geneCE[h])


what=matrix(nrow=1,ncol=15)




for (h in 1:nreps22) {
	
Emb6hr1 =subset(Emb6hr,V4>=startCE[h] & V5<=stopCE[h] )
EarlyAdultFemale1=subset(EarlyAdultFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyAdultMale1 =subset(EarlyAdultMale,V4>=startCE[h] & V5<=stopCE[h])
EarlyLastLarval1 =subset(EarlyLastLarval,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalFemale1 =subset(EarlyPupalFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalMale1=subset(EarlyPupalMale,V4>=startCE[h] & V5<=stopCE[h])
LateLastLarval1 =subset(LateLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MidLastLarval1 =subset(MidLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MIdPupalFemale1 =rep(NA,length(Emb6hr1[,1]))
MidPupalMale1 =subset(MidPupalMale,V4>=startCE[h] & V5<=stopCE[h])


what1=cbind(as.numeric(paste(Emb6hr1$V6)),as.numeric(paste(EarlyAdultFemale1$V6)),as.numeric(paste(EarlyAdultMale1$V6)),as.numeric(paste(EarlyLastLarval1$V6)),as.numeric(paste(EarlyPupalFemale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(LateLastLarval1$V6)),as.numeric(paste(MidLastLarval1$V6)),as.numeric(paste(MIdPupalFemale1)),as.numeric(paste(MidPupalMale1$V6)), startCE[h], stopCE[h], exonCE[h], geneCE[h])

print(h)

if (length(what1)>5){
what=rbind(what,what1)
}
}

whatPos=what


#####

Ch2pos=subset(Ch2, V7=="-")

PT=subset(Ch2pos, V3=="processed_transcript")
CE=subset(Ch2pos, V3=="coding_exon")


nreps=length(PT[,1])
start= numeric (nreps)
stop= numeric (nreps)
gene= numeric (nreps)

nreps2=length(CE[,1])
startCE= numeric (nreps2)
stopCE= numeric (nreps2)
exonCE= numeric (nreps2)
geneCE=numeric (nreps2)

for (i in 1:nreps) {
PT1=PT[i,]
start[i]=PT1$V4
stop[i]=PT1$V5	
gene[i]=i
}


k=1
num=0
i=1

for (i in 1:nreps2) {
CE1=subset(CE,V4>=start[i] & V5<=stop[i])
geneCCE=gene[i]
nreps3=length(CE1[,1])

for (k in 1:nreps3) {
num=num+1
startCE[num]=CE1$V4[k]
stopCE[num]=CE1$V5[k]
exonCE[num]=k
geneCE[num]=geneCCE
print(num)
}
}

xxx=cbind(startCE, stopCE, exonCE, geneCE)



Emb6hr=subset(Ch2, V2=="6hr_embryonic" & V3=="region")
EarlyAdultFemale =subset(Ch2, V2=="EarlyAdultFemale" & V3=="region")
EarlyAdultMale =subset(Ch2, V2=="EarlyAdultMale" & V3=="region")
EarlyLastLarval =subset(Ch2, V2=="EarlyLastLarval" & V3=="region")
EarlyPupalFemale =subset(Ch2, V2=="EarlyPupalFemale" & V3=="region")
EarlyPupalMale =subset(Ch2, V2=="EarlyPupalMale"& V3=="region")
LateLastLarval =subset(Ch2, V2=="LateLastLarval"& V3=="region")
MidLastLarval =subset(Ch2, V2=="MidLastLarval"& V3=="region")
MIdPupalFemale =rep(NA,length(Emb6hr1[,1]))
MidPupalMale =subset(Ch2, V2=="MidPupalMale"& V3=="region")



nreps22=length(xxx[,1])

h=1
	
Emb6hr1 =subset(Emb6hr,V4>=startCE[h] & V5<=stopCE[h])
EarlyAdultFemale1=subset(EarlyAdultFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyAdultMale1 =subset(EarlyAdultMale,V4>=startCE[h] & V5<=stopCE[h])
EarlyLastLarval1 =subset(EarlyLastLarval,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalFemale1 =subset(EarlyPupalFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalMale1=subset(EarlyPupalMale,V4>=startCE[h] & V5<=stopCE[h])
LateLastLarval1 =subset(LateLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MidLastLarval1 =subset(MidLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MIdPupalFemale1 =rep(NA,length(Emb6hr1[,1]))
MidPupalMale1 =subset(MidPupalMale,V4>=startCE[h] & V5<=stopCE[h])


#what=cbind(as.numeric(paste(Emb6hr1$V6)),as.numeric(paste(EarlyAdultFemale1$V6)),as.numeric(paste(EarlyAdultMale1$V6)),as.numeric(paste(EarlyLastLarval1$V6)),as.numeric(paste(EarlyPupalFemale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(LateLastLarval1$V6)),as.numeric(paste(MidLastLarval1$V6)),as.numeric(paste(MIdPupalFemale1$V6)),as.numeric(paste(MidPupalMale1$V6)), startCE[h], stopCE[h], exonCE[h], geneCE[h])




what=matrix(nrow=1,ncol=15)


for (h in 1:nreps22) {
	
Emb6hr1 =subset(Emb6hr,V4>=startCE[h] & V5<=stopCE[h] )
EarlyAdultFemale1=subset(EarlyAdultFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyAdultMale1 =subset(EarlyAdultMale,V4>=startCE[h] & V5<=stopCE[h])
EarlyLastLarval1 =subset(EarlyLastLarval,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalFemale1 =subset(EarlyPupalFemale,V4>=startCE[h] & V5<=stopCE[h])
EarlyPupalMale1=subset(EarlyPupalMale,V4>=startCE[h] & V5<=stopCE[h])
LateLastLarval1 =subset(LateLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MidLastLarval1 =subset(MidLastLarval,V4>=startCE[h] & V5<=stopCE[h])
MIdPupalFemale1 =rep(NA,length(Emb6hr1[,1]))
MidPupalMale1 =subset(MidPupalMale,V4>=startCE[h] & V5<=stopCE[h])


what1=cbind(as.numeric(paste(Emb6hr1$V6)),as.numeric(paste(EarlyAdultFemale1$V6)),as.numeric(paste(EarlyAdultMale1$V6)),as.numeric(paste(EarlyLastLarval1$V6)),as.numeric(paste(EarlyPupalFemale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(EarlyPupalMale1$V6)),as.numeric(paste(LateLastLarval1$V6)),as.numeric(paste(MidLastLarval1$V6)),as.numeric(paste(MIdPupalFemale1)),as.numeric(paste(MidPupalMale1$V6)), startCE[h], stopCE[h], exonCE[h], geneCE[h])

print(h)

if (length(what1)>5){
what=rbind(what,what1)
}
}


whatNeg=what
whatNeg[,15]=whatNeg[,15]+.5


what2=rbind(whatPos,whatNeg)


lineplot.CI(what2[,15],(what2[,1]),type="p")


data=as.data.frame((what2))

#####################   CHANGE ME    #####################
#####################   CHANGE ME    #####################

data$Chr="Ch9"

#####################   CHANGE ME    #####################
#####################   CHANGE ME    #####################


data$x=paste(data$Chr,data$V12,data$V13,data$V15, sep="_")

EarlyAdultFemale1yBV=data$V1-data$V2
EarlyAdultFemale1xBV=data$x
EarlyAdultFemale1Data=cbind.data.frame(EarlyAdultFemale1yBV, EarlyAdultFemale1xBV,c("EarlyAdultFemale1"))
names(EarlyAdultFemale1Data)=c("y","x","fact")

EarlyAdultMale1yBV=data$V1-data$V3
EarlyAdultMale1xBV=data$x
EarlyAdultMale1Data=cbind.data.frame(EarlyAdultMale1yBV, EarlyAdultMale1xBV,"EarlyAdultMale1")
names(EarlyAdultMale1Data)=c("y","x","fact")

EarlyLastLarval1yBV=data$V1-data$V4
EarlyLastLarval1xBV=data$x
EarlyLastLarval1Data=cbind.data.frame(EarlyLastLarval1yBV, EarlyLastLarval1xBV,"EarlyLastLarval1")
names(EarlyLastLarval1Data)=c("y","x","fact")

EarlyPupalFemale1yBV=data$V1-data$V5
EarlyPupalFemale1xBV=data$x
EarlyPupalFemale1Data=cbind.data.frame(EarlyPupalFemale1yBV, EarlyPupalFemale1xBV,"EarlyPupalFemale1")
names(EarlyPupalFemale1Data)=c("y","x","fact")

EarlyPupalMale1yBV=data$V1-data$V6
EarlyPupalMale1xBV=data$x
EarlyPupalMale1Data=cbind.data.frame(EarlyPupalMale1yBV, EarlyPupalMale1xBV,"EarlyPupalMale1")
names(EarlyPupalMale1Data)=c("y","x","fact")

LateLastLarval1yBV=data$V1-data$V8
LateLastLarval1xBV=data$x
LateLastLarval1Data=cbind.data.frame(LateLastLarval1yBV, LateLastLarval1xBV,"LateLastLarval1")
names(LateLastLarval1Data)=c("y","x","fact")

MidLastLarval1yBV=data$V1-data$V9
MidLastLarval1xBV=data$x
MidLastLarval1Data=cbind.data.frame(MidLastLarval1yBV, MidLastLarval1xBV,"MidLastLarval1")
names(MidLastLarval1Data)=c("y","x","fact")

MIdPupalFemale1yBV=data$V1-data$V10
MIdPupalFemale1xBV=data$x
MIdPupalFemale1Data=cbind.data.frame(MIdPupalFemale1yBV, MIdPupalFemale1xBV,"MIdPupalFemale1")
names(MIdPupalFemale1Data)=c("y","x","fact")

MidPupalMale1yBV=data$V1-data$V11
MidPupalMale1xBV=data$x
MidPupalMale1Data=cbind.data.frame(MidPupalMale1yBV, MidPupalMale1xBV,"MidPupalMale1")
names(MidPupalMale1Data)=c("y","x","fact")

grandMData=rbind(EarlyAdultFemale1Data, EarlyAdultMale1Data, EarlyLastLarval1Data, EarlyPupalFemale1Data, EarlyPupalMale1Data, LateLastLarval1Data, MidLastLarval1Data, MIdPupalFemale1Data, MidPupalMale1Data)


#z=lm(y~x+fact,data=grandMData)

z1=lme(y~ fact,random=~1|x,data=grandMData,na.action = "na.omit")

zz=(coef(z1))
zzq=quantile(zz[,1],.99)
zzqs=subset(zz,zz[,1]>=zzq)
zzqs
max(zzqs)
head(grandMData)
tail(grandMData)

write.csv(grandMData, file="ChLG9_20000000_21459655.csv")






#################################

x1=read.csv("ChLG2_1_4000000.csv")
x2=read.csv("ChLG2_4000000_8000000.csv")
x3=read.csv("ChLG2_8000000_12000000.csv")
x4=read.csv("ChLG2_12000000_16000000.csv")
x5=read.csv("ChLG2_16000000_20218415.csv")
x6=read.csv("ChLG3_1_4000000.csv")
x7=read.csv("ChLG3_4000000_8000000.csv")
x8=read.csv("ChLG3_8000000_12000000.csv")
x9=read.csv("ChLG3_12000000_16000000.csv")
x10=read.csv("ChLG3_16000000_20000000.csv")
x11=read.csv("ChLG3_20000000_24000000.csv")
x12=read.csv("ChLG3_24000000_28000000.csv")
x13=read.csv("ChLG3_28000000_32000000.csv")
x14=read.csv("ChLG3_32000000_36000000.csv")
x15=read.csv("ChLG3_36000000_38791480.csv")
x16=read.csv("ChLG4_1_4000000.csv")
x17=read.csv("ChLG4_4000000_8000000.csv")
x18=read.csv("ChLG4_8000000_12000000.csv")
x19=read.csv("ChLG4_12000000_13894384.csv")
x20=read.csv("ChLG5_1_4000000.csv")
x21=read.csv("ChLG5_4000000_8000000.csv")
x22=read.csv("ChLG5_8000000_12000000.csv")
x23=read.csv("ChLG5_12000000_16000000.csv")
x24=read.csv("ChLG5_16000000_19135781.csv")
x25=read.csv("ChLG6_1_4000000.csv")
x26=read.csv("ChLG6_4000000_8000000.csv")
x27=read.csv("ChLG6_8000000_12000000.csv")
x28=read.csv("ChLG6_12000000_13176827.csv")
x29=read.csv("ChLG7_1_4000000.csv")
x30=read.csv("ChLG7_4000000_8000000.csv")
x31=read.csv("ChLG7_8000000_12000000.csv")
x32=read.csv("ChLG7_12000000_16000000.csv")
x33=read.csv("ChLG7_16000000_20532854.csv")
x34=read.csv("ChLG8_1_4000000.csv")
x35=read.csv("ChLG8_4000000_8000000.csv")
x36=read.csv("ChLG8_8000000_12000000.csv")
x37=read.csv("ChLG8_12000000_16000000.csv")
x38=read.csv("ChLG8_16000000_18021898.csv")
x39=read.csv("ChLG9_1_4000000.csv")
x40=read.csv("ChLG9_4000000_8000000.csv")
x41=read.csv("ChLG9_8000000_12000000.csv")
x42=read.csv("ChLG9_12000000_16000000.csv")
x43=read.csv("ChLG9_16000000_20000000.csv")
x44=read.csv("ChLG9_20000000_21459655.csv")
x45=read.csv("ChLG10_1_4000000.csv")
x46=read.csv("ChLG10_4000000_8000000.csv")
x47=read.csv("ChLG10_8000000_11386040.csv")
x48=read.csv("ChLGX_1_4000000.csv")
x49=read.csv("ChLGX_4000000_8000000.csv")
x50=read.csv("ChLGX_8000000_10877635.csv")

BAM=rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,x41,x42,x43,x44,x45,x46,x47,x48,x49,x50)


z1=lme(y~ fact,random=~1|x,data=BAM,na.action = "na.omit")

zz=(coef(z1))
zzq=quantile(zz[,1],.995)
zzqs=subset(zz,zz[,1]>=zzq)
zzqs
max(zzqs)
head(BAM)
tail(BAM)


write.csv(zzqs, file="BVs.csv")


