#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)


#Make sure that there is no open graphs
graphics.off()


#Read the simulations results
cluster200<-read.table("../additional_simulations/clusters_200.dat",header=TRUE)
cluster200$target<-"n[target]==200"

cluster104<-read.table("../additional_simulations/clusters_10_4.dat",header=TRUE)
cluster104$target<-"n[target]==10^4"


#Definition of characteristics of the figure
font = 24
w=1
wid_line=0.7
wid_comet=4


#Produce the distribution of clusters plots
p1<-ggplot(cluster200,aes(x=cluster, y=amount))+geom_bar(stat = "identity")+theme_bw(base_size=font)+xlab("Cluster size")+ylab("Number of clusters per cell")+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9),labels=c("1","2","3","4","5","6","7","8","9")) + coord_cartesian(xlim=c(1,7))
p3<-ggplot(cluster104,aes(x=cluster, y=amount))+geom_bar(stat = "identity")+theme_bw(base_size=font)+xlab("Cluster size")+ylab("Number of clusters per cell")+scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9),labels=c("1","2","3","4","5","6","7","8","9")) + coord_cartesian(xlim=c(1,7))


#Generate comet plots
comet<-matrix(data = scan("../data/abrams_comet_binning_40.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
br=seq(0,40, by=4)

dist_sim_200 <- c(0,cluster200$population)
dist_sim_104 <- c(0,cluster104$population)
dist_exp <- comet[1,] #At 15 min
y <- c()
x_sim_200 <- c()
x_sim_104 <- c()
x_exp <- c()

for (i in 1:(length(dist_exp))) {
  y<- c(y, br[i], br[i+1])
  
  x_sim_200 <- c(x_sim_200, dist_sim_200[i], dist_sim_200[i])
  x_sim_104 <- c(x_sim_104, dist_sim_104[i], dist_sim_104[i])
  x_exp <- c(x_exp, dist_exp[i], dist_exp[i])
}

y<-c(y, rev(y), 0)
x_sim_200 <- c(-x_sim_200, rev(x_sim_200), -x_sim_200[1])
x_sim_104 <- c(-x_sim_104, rev(x_sim_104), -x_sim_104[1])
x_exp <- c(2-x_exp, 2+rev(x_exp), 2-x_exp[1])


t<-c(rep(0, length(x_sim_200)),rep(2, length(x_sim_200)))
amount_200<-c(x_sim_200,x_exp)
amount_104<-c(x_sim_104,x_exp)
tail<-rep(y,2)

comet_200 <- data.frame(t,amount_200,tail)
names(comet_200) <- c('t', 'amount', 'tail')
comet_200$target<-"n[target]==200"

comet_104 <- data.frame(t,amount_104,tail)
names(comet_104) <- c('t', 'amount', 'tail')
comet_104$target<-"n[target]==10^4"


#Produce the synthetic comet plot 
p2<-ggplot(comet_200, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line)+xlab("Source")+ylab("Cell franction per population")+theme_bw(base_size=font)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(rows =vars(target),labeller = label_parsed)+scale_y_continuous(breaks=seq(2,38, by=4), labels = c(expression(L[0]),expression(L[1]),expression(L[2]),expression(L[3]),expression(L[4]),expression(L[5]),expression(L[6]),expression(L[7]),expression(L[8]),expression(L[9])))+scale_x_continuous(breaks=c(0,2),labels=c("Simulation","Experiment"))
p4<-ggplot(comet_104, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line)+xlab("Source")+ylab("Cell franction per population")+theme_bw(base_size=font)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(rows =vars(target),labeller = label_parsed)+scale_y_continuous(breaks=seq(2,38, by=4), labels = c(expression(L[0]),expression(L[1]),expression(L[2]),expression(L[3]),expression(L[4]),expression(L[5]),expression(L[6]),expression(L[7]),expression(L[8]),expression(L[9])))+scale_x_continuous(breaks=c(0,2),labels=c("Simulation","Experiment"))


#Save the figure
pdf("Fig_8.pdf",width=16,height=12)
ggarrange(p1+theme(axis.text.x=element_blank())+xlab(""), p2+theme(axis.text.x=element_blank())+xlab(""), p3, p4, labels = c("A", "B","C","D"), font.label = list(size = font), ncol = 2, nrow = 2)
graphics.off()