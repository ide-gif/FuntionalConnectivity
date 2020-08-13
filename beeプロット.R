#.libPaths("C:/Users/nryu-/AppData/Roaming/R")

#パッケージの読み込み
library("ggbeeswarm")
library("ggplot2")
###データ例の作成#####
data.frame1=read.delim("clipboard")

#######

#プロット雛形を作成
gg <- ggplot()
gg=gg+theme_classic()
#Beeswarmプロット:geom_quasirandomコマンド
gg=gg + geom_quasirandom(data = data.frame1, 
                         aes(x = factor(DIV), 
                             y = value,col=type),
                         stat = "identity",
                         dodge.width = 0.7, cex = 2)
gg=gg+xlab("培養日数")
gg=gg+ylab("平均に対する発火頻度")
gg=gg+theme(axis.title = element_text(size=20),
            axis.text = element_text(size=14,
                                     colour = "black"),
            legend.text = element_text(size = 15))
#gg=gg+geom_hline(yintercept=1,colour="red",size=2)
gg

#400*400