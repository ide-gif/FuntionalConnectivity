rm(list = ls())
setwd("C:/Ide_AQUACOSMOS/191228AAV-2/area2/DIV13")
x<-read.csv("Results.csv")
Size<-dim(x)
Time=Size[1]
Back=Size[2]
ROI=Back-2 #TimeとBack分を引く#
i=2
Int=array(dim=c(Time,ROI))
#SubtractBackground#
while(i<Back){
  j=1
  while(j<=Time){
    Int[j,i-1]=x[j,i]-x[j,Back]
    j=j+1
  }
  i=i+1
}


#F/F0#
F0=Int[1,]
i=1
while (i<=ROI) {
  j=1
  while(j<=Time){
    Int[j,i]=Int[j,i]/F0[i]
    j=j+1
  }
  i=i+1
}

#微分の平均とSEMをとる#
#微分#
D=array(dim=c(Time,ROI))
i=1
j=1
while (i<=ROI) {
  j=1
  while (j<Time) {
    D[j,i]=((Int[j+1,i]-Int[j,i])/(x[j+1,1]-x[j,1]))
    j=j+1
  }
  i=i+1
}

i=1
while (i<=ROI) {
  D[Time,i]=D[Time-1,i]
  i=i+1
}

#微分の平均#
Dave=array(dim=c(ROI))
i=1
while (i<=ROI) {
  Dave[i]=mean(D[,i])
  i=i+1
}

#微分の標準偏差#
DSEM=array(dim=c(ROI))
i=1
while (i<=ROI) {
  DSEM[i]=sd(D[,i])
  i=i+1
}

#閾値の設定#
Dshre=array(dim=c(ROI))
i=1
while (i<=ROI) {
  Dshre[i]=Dave[i]+DSEM[i]*3
  i=i+1
}

#ラスター1,0で表す#
ras=array(dim=c(ROI,Time))
i=1
j=1
while (i<=ROI) {
  j=1
  while (j<=Time) {
    if(D[j,i]>Dshre[i]){
      ras[i,j]=1
      j=j+1
    }
    else{
      ras[i,j]=0
      j=j+1
    }
  }
  i=i+1
}

#ラスターの1が続いているところを除去する#
i=1
j=Time
while (i<=ROI) {
  j=Time
  while (j>1) {
    if((ras[i,j]==1)&&(ras[i,j-1]==1)){
      ras[i,j]=0
      j=j-1
    }
    else{
      j=j-1
    }
  }
  i=i+1
}

#1細胞の総発火数
fire=array(dim=c(ROI))
i=1
j=1
fire[]=0
while (i<=ROI) {
  j=1
  while (j<=Time) {
    fire[i]=fire[i]+ras[i,j]
    j=j+1
  }
  i=i+1
}

#T1min当たりの発火数にする
fire[]=fire[]/17