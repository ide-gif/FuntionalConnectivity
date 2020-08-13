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
ymax=0.1
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

f3="出力スケールフリー.png"
png(f3,width=350,height = 400)
#グラフの調整
par(mar=c(5,6,1,1),
    mgp=c(3.8,0.7,0),
    lwd=3,
    tcl=0.4,
    ps=25,
    las=1,
    bty="o",
    xaxs="r",
    yaxs="r"
)
par(mfrow=c(1,1))
plot(oc*100,op,log="xy",xlim=c(xmin,xmax),ylim=c(ymin,ymax),ann=F,cex=1.5,col="black",pch=16,cex.axis=0.8)
par(new=T)
plot(z0*100,z1,type="l",log="xy",xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab = "全細胞に対する出力の割合（％）",ylab="存在確率",
     lwd=4,col="black",cex.axis=0.8
)
dev.off()

#polygon(polyx,polyy,col="#FF00007F",lwd=0.01^100)

f4="入力スケールフリー.png"
png(f4,width=350,height = 400)
#グラフの調整
par(mar=c(5,6,1,1),
    mgp=c(3.8,0.7,0),
    lwd=3,
    tcl=0.4,
    ps=25,
    las=1,
    bty="o",
    xaxs="r",
    yaxs="r"
)
par(mfrow=c(1,1))
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

