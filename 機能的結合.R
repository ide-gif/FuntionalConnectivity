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
