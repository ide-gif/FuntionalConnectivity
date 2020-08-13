
###データ例の作成#####
#data.frame2=read.delim("clipboard")

#######
#Beeswarmプロット:geom_quasirandomコマンド
gg=gg + geom_quasirandom(data = data.frame2, 
                             aes(x = factor(DIV), 
                                 y = value,col=type),
                             stat = "identity",
                             dodge.width = 0.7, cex = 2)
gg=gg+xlab("培養日数")
gg=gg+ylab("平均結合数に対する結合数の割合")


gg