#Open the required libraries
library(reshape2)
library(ggplot2)


#Characteristics of the figure
size=30
wid_line=1


#Read the simulation results
data<-read.table("../additional_simulations/survival_after_dr_discriminators.txt",header=TRUE)
S<-c(data$S_001,data$S_01,data$S_2,data$S_20)
dr <-c(rep(0.01,length(data$S_001)),rep(0.1,length(data$S_001)),rep(2,length(data$S_001)),rep(20,length(data$S_001)))
D<-rep(data$dose,4)
survival<-data.frame(D,dr,S)
names(survival) <- c('dose', 'dose_rate', 'S')


#Produce and save the figure
pdf("Fig_5.pdf",width=12,height=8)
print(
  ggplot(survival, aes(x=dose, y=S, color=factor(dose_rate), group=dose_rate))+geom_line(size=wid_line)+
    scale_color_manual(values=c("orange", "blue", "red", "black"), labels=c("0.01 Gy/min", "0.1 Gy/min", "2 Gy/min", "20 Gy/min"))+
    xlab("Dose / Gy")+ylab("Surviving fraction of cells")+
    scale_y_log10(breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1))+
    theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position=c(0.2,0.2),legend.title=element_blank(),legend.text = element_text(size = size),legend.box.background=element_rect(linewidth = 2),legend.margin=margin(0.5,3,3,3, unit='mm'))
)

graphics.off()