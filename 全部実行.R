rm(list = ls())
setwd("C:/Ide_AQUACOSMOS/200117/DIV8-fluo4-2")
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

i=1
FC=array(dim=c(ROI,ROI)) #機能的結合をもつとき値を入れる#
FC[]=0
while(i<=ROI){
  k=1
  while (k<=ROI) {
    j=5
    distribution=c() #細胞Aの発火に対する細胞Bの発火の分布#
    while (j<=Time-4) {
      if(ras[i,j]==1){
        if(ras[k,j-4]==1){
          distribution=c(distribution,-4)
        }
        if(ras[k,j-3]==1){
          distribution=c(distribution,-3)
        }
        if(ras[k,j-2]==1){
          distribution=c(distribution,-2)
        }
        if(ras[k,j-1]==1){
          distribution=c(distribution,-1)
        }
        if(ras[k,j]==1){
          distribution=c(distribution,0)
        }
        if(ras[k,j+1]==1){
          distribution=c(distribution,1)
        }
        if(ras[k,j+2]==1){
          distribution=c(distribution,2)
        }
        if(ras[k,j+3]==1){
          distribution=c(distribution,3)
        }
        if(ras[k,j+4]==1){
          distribution=c(distribution,4)
        }
      }
      j=j+1
    }
    
    #検定#
    ks=tryCatch({
      ks.test(x=distribution,y="pnorm",mean=0,sd=sd(distribution),exact=T)$p.value
    },
    error=function(err){
      return(1)
    })
    
    #t検定は最低3以上のベクトルが必要
    t=tryCatch({
      t.test(distribution,mu=0)$p.value
    },
    error=function(err){
      return(1)
    },
    warning=function(warn){
      return(1)
    }) #tはt検定のp値を算出している。エラー分が出るようなときはt=1を返すようにしている#
    ks=ifelse(is.na(ks),1,ks)
    t=ifelse(is.na(t),1,t)
    
    #FCに値を入れる
    if((ks<=0.05)&&(t<=0.05)){
      if(mean(distribution)<0){
        #kの後にiが発火する
        FC[k,i]=1
      }
      else if(mean(distribution)>=0){
        #iの後にkが発火する
        FC[i,k]=1
      }
    }
    #FC[i,j]=1はi→jを意味する
    k=k+1
  }
  i=i+1
}


#グラフの調整
par(mar=c(0.5,5,1,0.5),
    mgp=c(3,0.5,0),
    bty="o",
    lwd=3,
    tcl=0.4,
    ps=18,
    las=1,
    mfrow=c(1,1)
)
sink("発火数.txt")
#発火数
print("発火数")
fire
#発火多い順
print("発火ランキング")
order(-fire)
#発火の平均
print("平均発火数")
mean(fire)
sink()
#hist(fire,xlab="発火頻度 (/min)",ylab="細胞数の割合",freq = F,breaks=6,main="")
#横400縦250
f1="発火した細胞の割合.png"
png(f1,width =400,height = 250 )
par(mar=c(0.5,5,1,0.5),
    mgp=c(3,0.5,0),
    bty="o",
    lwd=3,
    tcl=0.4,
    ps=18,
    las=1,
    mfrow=c(1,1)
)
barplot(colSums(ras)/ROI,ylim=c(0,1),xlab = "",ylab = "発火した細胞の割合",cex.axis=0.8)
dev.off()

#20%以上の細胞が同時発火したタイミングが何回あるか
burst=array(dim=c(Time))
burst=colSums(ras)/ROI
nburst=0
i=1
while(i<=Time){
  if(burst[i]>0.2){
    nburst=nburst+1
  }
  i=i+1
}
nburst

#グラフの調整
par(mar=c(5,5,1,1),
    mgp=c(3,0.7,0),
    lwd=3,
    tcl=0.4,
    ps=18,
    las=1,
    bty="o",
    xaxs="r",
    yaxs="r"
)

#outputしている数
rS=rowSums(FC)
rS
#input
cS=colSums(FC)
cS
order(-rS)
order(-cS)
##ROIごとにおける出力数＋入力数の数をSSとする
SS=rS+cS
SS
consum=c()
#SSの0部分を削除する
i=1
while (i<=ROI) {
  if(SS[i]!=0){
    consum=c(consum,SS[i])
  }
  i=i+1
}
consum

numberofactive=length(consum)


#横800縦500のサイズが良い
#発火数と結合数の関係
f2="発火数と結合数.png"
png(f2,width = 800,height = 500)
#グラフの調整
par(mar=c(5,5,1,1),
    mgp=c(3,0.7,0),
    lwd=3,
    tcl=0.4,
    ps=18,
    las=1,
    bty="o",
    xaxs="r",
    yaxs="r"
)
par(mfrow=c(1,2))
plot(fire,rS,xlab = "発火数（/分）",ylab="出力数",pch=16,cex=1.5,cex.axis=0.8)
plot(fire,cS,xlab = "発火数（/分）",ylab="入力数",pch=16,cex=1.5,cex.axis=0.8)
par(mfrow=c(1,1))
dev.off()
#相関係数の算出
cor.test(fire,rS)
cor.test(fire,cS)

#fitting curveを描く
i=1
#結合の存在確率を表すベクトル（縦軸）
outprob=c()
inprob=c()
#結合を表すベクトル（横軸）
outcon=c(1:max(rS))
incon=c(1:max(cS))
while (i<=max(rS)) {
  outprob=c(outprob,length(rS[i==rS]))
  i=i+1
}
i=1
while (i<=max(cS)) {
  inprob=c(inprob,length(cS[i==cS]))
  i=i+1
}


#実際に使うベクトル
op=c()
ip=c()
oc=c()
ic=c()

#outprobの0部分を削除する
i=1
while (i<=max(rS)) {
  if(outprob[i]!=0){
    op=c(op,outprob[i])
    oc=c(oc,i)
  }
  i=i+1
}
#inprobの0部分を削除する
i=1
while (i<=max(cS)) {
  if(inprob[i]!=0){
    ip=c(ip,inprob[i])
    ic=c(ic,i)
  }
  i=i+1
}
op
oc
ip
ic
#結合数iの細胞/結合をもつ全細胞
op=op/numberofactive
ip=ip/numberofactive
#全細胞の何％と結合しているか
oc=oc/numberofactive
ic=ic/numberofactive

#軸の調整
xmin=2
xmax=100
ymin=min(min(ip),min(op))
ymax=0.2
outkinji=nls(op~a*(oc^(b)),start=c(a=0.1,b=(-1)))
outkinji
a1=coef(outkinji)[1]
b1=coef(outkinji)[2]
z0=seq(0.01,1,by=0.05)
z1=a1*(z0^(b1))

#high connectivityの領域設定
polyx=c(40,40,115,115)
polyy=c(0.0000001,10,10,0.0000001)

#input
inkinji=nls(ip~a*(ic^(b)),start=c(a=0.1,b=(-1)))
inkinji
a2=coef(inkinji)[1]
b2=coef(inkinji)[2]
z2=seq(0.01,1,by=0.05)
z3=a2*(z2^(b2))

f3="スケールフリー.png"
png(f3,width=800,height = 500)
#グラフの調整
par(mar=c(5,5,1,1),
    mgp=c(3,0.7,0),
    lwd=3,
    tcl=0.4,
    ps=18,
    las=1,
    bty="o",
    xaxs="r",
    yaxs="r"
)
par(mfrow=c(1,2))
plot(oc*100,op,log="xy",xlim=c(xmin,xmax),ylim=c(ymin,ymax),ann=F,cex=1.5,col="black",pch=16,cex.axis=0.8)
par(new=T)
plot(z0*100,z1,type="l",log="xy",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab = "全細胞に対する出力の割合（％）",ylab="存在確率",
     lwd=4,col="black",cex.axis=0.8
)
#polygon(polyx,polyy,col="#FF00007F",lwd=0.01^100)
par(new=F)
plot(ic*100,ip,log="xy",xlim=c(xmin,xmax),ylim=c(ymin,ymax),ann=F,cex=1.5,col="black",pch=16,cex.axis=0.8)
par(new=T)
plot(z2*100,z3,type="l",log="xy",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab = "全細胞に対する入力の割合（％）",ylab="存在確率",
     lwd=4,col="black",cex.axis=0.8)
#polygon(polyx,polyy,col="#FF00007F",lwd=0.01^100)
par(mfrow=c(1,1))
dev.off()
highshre=0.4*numberofactive
sink("機能的結合.txt")
#出力
print("出力")
rS
print("出力ランキング")
order(-rS)
#入力
print("入力")
cS
print("入力ランキング")
order(-cS)
#発火数と結合数の相関係数
print("発火数と結合数の相関係数")
print("出力")
cor.test(fire,rS)
print("入力")
cor.test(fire,cS)
#スケールフリーの近似曲線
print("スケールフリーの近似曲線")
print("出力")
outkinji
print("入力")
inkinji
#平均結合数
print("平均結合数")
sum(cS)/ROI
#0.4以上でhigh connectivity
print("0.4以上でhigh connectivity")
print("高出力ニューロン")
(1:length(rS))[rS>highshre]
print("高入力ニューロン")
(1:length(cS))[cS>highshre]
sink()

#結合をもつ細胞の出力
sink("結合の方向性.txt")
i=1
while (i<=ROI) {
  j=1
  while (j<=ROI) {
    if(FC[i,j]==1){
      cat(i,"→",j,"\n")
    }
    j=j+1
  }
  i=i+1
}
sink()