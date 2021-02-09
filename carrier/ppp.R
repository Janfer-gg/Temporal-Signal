library(ggplot2)

df<-data.frame(x=1:9232,y=1)
p<-ggplot(df,aes(x=x,y=y))+ geom_path(color = "grey", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)

df1<-data.frame(x=3:229,y=1)
df2<-data.frame(x=230:410,y=1)
df3<-data.frame(x=457:582,y=1)
df4<-data.frame(x=1075:1308,y=1)
df5<-data.frame(x=1803:1920,y=1)
df6<-data.frame(x=1927:2175,y=1)
df7<-data.frame(x=2176:3097,y=1)
df8<-data.frame(x=3098:3173,y=1)
df9<-data.frame(x=3186:3690,y=1)
df10<-data.frame(x=3703:5082,y=1)
df11<-data.frame(x=5093:5681,y=1)
df12<-data.frame(x=5753:5986,y=1)
df13<-data.frame(x=6058:6179,y=1)
df14<-data.frame(x=7147:8007,y=1)
df15<-data.frame(x=8178:8766,y=1)

p<-p+geom_line(data=df1,aes(x=x,y=y),color="green",size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df2,aes(x=x,y=y),color="red",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df3,aes(x=x,y=y),color="yellow",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df4,aes(x=x,y=y),color="blue",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df5,aes(x=x,y=y),color="orange",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df6,aes(x=x,y=y),color="purple",alpha=.5,size=4,arrow=arrow(length = unit(0.20, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df7,aes(x=x,y=y),color="pink",alpha=.5,size=4,arrow=arrow(length = unit(0.20, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df8,aes(x=x,y=y),color="green",alpha=.5,size=4,arrow=arrow(length = unit(0.20, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df9,aes(x=x,y=y),color="red",alpha=.5,size=4,arrow=arrow(length = unit(0.20, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df10,aes(x=x,y=y),color="yellow",alpha=.5,size=4,arrow=arrow(length = unit(0.20, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df11,aes(x=x,y=y),color="blue",alpha=.5,size=4,arrow=arrow(length = unit(0.20, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df12,aes(x=x,y=y),color="orange",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df13,aes(x=x,y=y),color="purple",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df14,aes(x=x,y=y),color="pink",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
  geom_line(data=df15,aes(x=x,y=y),color="blue2",alpha=.5,size=4,arrow=arrow(length = unit(0.30, "cm"),ends = "last",type = "closed"))+
coord_polar(theta = "x")

png("C://Users//41518//Desktop//111.png",width = 480*3,height = 480*3,res = 3)
print(p)
dev.off()


p+annotate("segment",x=1,xend = 2000,y=1,yend=1,arrow=arrow(angle=40,length = unit(0.4, "cm"),ends = "first",type = "closed"),size=2,color="black")

# p+geom_line(data=df1,aes(x=x,y=y),color="black",alpha=.5,size=5)+geom_segment(x=3,xend = 2000,yend=1,arrow=arrow(length = unit(0.80, "cm"),ends = "last",type = "closed"),color="black")


