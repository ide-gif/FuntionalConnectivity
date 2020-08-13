distribution=c() #細胞Aの発火に対する細胞Bの発火の分布#
i=2
k=29
j=5
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


ks=ks.test(x=distribution,y="pnorm",mean=0,sd=sd(distribution),exact=T)$p.value
t=t.test(distribution,mu=0)$p.value
#tはt検定のp値を算出している。エラー分が出るようなときはt=1を返すようにしている#
ks
t
distribution
mean(distribution)
hist(distribution)