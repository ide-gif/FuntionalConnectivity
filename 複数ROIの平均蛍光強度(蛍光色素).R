rm(list = ls())
name="‚ "
x<-read.delim("clipboard")
Size<-dim(x)
Time=Size[1]
Back=Size[2]
ROI=Back-2 #Time‚ÆBack•ª‚ðˆø‚­#
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

SEM=array(dim=c(Time))
j=1
while(j<=Time){
  SEM[j]=sd(Int[j,])/sqrt(ROI)
  j=j+1
}
Av=rowMeans(Int)
Ti=x[,1]
library(ggplot2)
Y=data.frame(
  T=x[,1],Ave=rowMeans(Int),SEM
)
g=ggplot(Y,aes(T,Ave))
g=g+ggtitle(name)
g=g+geom_linerange(aes(ymin=Ave-SEM,ymax=Ave+SEM,colour="red"))
g=g+geom_line(size=1.5)
g=g+guides(colour=FALSE)
g=g+labs(x="Time (sec)",y="ƒ¢F/F0")
g=g+theme(
  axis.line = element_line("black"),
  axis.ticks = element_line("black"),
  axis.text.x=element_text(margin=unit(rep(8,4),"mm")),
  axis.text.y=element_text(margin=unit(rep(8,4),"mm")),
  axis.ticks.length=unit(-3,"mm"),
  axis.text = element_text(size=12,face=1),
  axis.title = element_text(size=14,face=1),
  panel.background = element_blank())
g=g+scale_x_continuous(expand = c(0, 0))
g
extension=".png"
filename=paste(name,extension,sep="")
#ggsave(file=filename)#