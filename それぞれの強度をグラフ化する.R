
i=38
g=i+8
par(mfrow=c(3,3))
while (i<=g) {
  plot(x[,1],Int[,i],type = "l",ann = F)
  i=i+1
  }
par(mfrow=c(1,1))

#plot(x[,1],Int[,i],type = "l",ann = F)
fire[i-9]*17