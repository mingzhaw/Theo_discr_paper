#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)


#Make sure that there is no open graphs
graphics.off()


#Definition of characteristics of the figure
w=1
size=24
wid_line=0.7
wid_clon=0.2
wid_comet=4
x11()


#Generation of synthetic comet from simulations of different alpha values
comet<-matrix(data = scan("../additional_simulations/comet_alpha.dat", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
br=seq(0,40, by=4)

dist_03 <- comet[1,]
dist_10 <- comet[2,]
dist_15 <- comet[3,]
dist_20 <- comet[4,]

y <- c()
x_03 <- c()
x_10 <- c()
x_15 <- c()
x_20 <- c()

for (i in 1:(length(dist_03))) {
  y<- c(y, br[i], br[i+1])
  
  x_03 <- c(x_03, dist_03[i], dist_03[i])
  x_10 <- c(x_10, dist_10[i], dist_10[i])
  x_15 <- c(x_15, dist_15[i], dist_15[i])
  x_20 <- c(x_20, dist_20[i], dist_20[i])
}

y<-c(y, rev(y), 0)
x_03 <- c(-x_03, rev(x_03), -x_03[1])
x_10 <- c(2-x_10, 2+rev(x_10), 2-x_10[1])
x_15 <- c(4-x_15, 4+rev(x_15), 4-x_15[1])
x_20 <- c(6-x_20, 6+rev(x_20), 6-x_20[1])


t<-c(rep(0, length(x_03)),rep(2, length(x_03)),rep(4, length(x_03)),rep(6, length(x_03)))
amount<-c(x_03,x_10,x_15,x_20)
tail<-rep(y,4)

comet_alpha <- data.frame(t,amount,tail)
names(comet_alpha) <- c('t', 'amount', 'tail')


#Read experimental comet data
data<-matrix(data = scan("../data/abrams_comet_binning_40.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
t <- c(0, 2, 4, 6)
time<-c()
amount<-c()
mean_tail<-rep(seq(2,38, by=4),length(t))
for (i in 1:(length(t))) {
  time<-c(time,rep(t[i],length(data[i,])))
  amount<-c(amount,data[1,]*w)
}
comet_real_data<- data.frame(time,amount,mean_tail)
names(comet_real_data) <- c('t', 'amount', 'tail')


#Save the figure
pdf("Fig_9.pdf",width=8,height=8)
ggplot(comet_alpha, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line, color="orange")+xlab(expression(alpha~"/"~Gy^-1))+ylab("Relative tail intensity / %")+geom_errorbar(data=comet_real_data,aes(xmin=t-amount,xmax=t+amount,y=tail),width=wid_comet,position=position_dodge(0.05),size=wid_line, color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+scale_y_continuous(breaks=c(0,10,20,30,40))+scale_x_continuous(breaks=c(0,2,4,6),labels=c("0.3","1.0","1.5","2.0"))
graphics.off()

