#Open required libraries
library(reshape2)
library(ggplot2)
library(ggpubr)

#Make sure that there is no open graphs
graphics.off()


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


#Function to generate violin plots from the simulations
comet_violin <- function(comet){
  br=seq(0,40, by=4)
  
  
  dist15 <- comet[1,]
  dist30 <- comet[2,]
  dist60 <- comet[3,]
  dist120 <- comet[4,]
  dist240 <- comet[5,]
  dist360 <- comet[6,]
  
  
  dist15 <- dist15*w
  dist30 <- dist30*w
  dist60 <- dist60*w
  dist120 <- dist120*w
  dist240 <- dist240*w
  dist360 <- dist360*w
  
  t <- c(0, 2, 4, 6, 8, 10)
  
  
  y <- c()
  x0 <- c()
  x15 <- c()
  x30 <- c()
  x60 <- c()
  x120 <- c()
  x240 <- c()
  x360 <- c()
  
  ymax=44
  
  for (i in 1:(length(dist15))) {
    y<- c(y, br[i], br[i+1])
    
    x15 <- c(x15, dist15[i], dist15[i])
    x30 <- c(x30, dist30[i], dist30[i])
    x60 <- c(x60, dist60[i], dist60[i])
    x120 <- c(x120, dist120[i], dist120[i])
    x240 <- c(x240, dist240[i], dist240[i])
    x360 <- c(x360, dist360[i], dist360[i])
  }
  
  y<-c(y, rev(y), 0)
  x15 <- c(-x15, rev(x15), -x15[1])
  x30 <- c(2-x30, 2+rev(x30), 2-x30[1])
  x60 <- c(4-x60, 4+rev(x60), 4-x60[1])
  x120 <- c(6-x120, 6+rev(x120), 6-x120[1])
  x240 <- c(8-x240, 8+rev(x240), 8-x240[1])
  x360 <- c(10-x360, 10+rev(x360), 10-x360[1])
  
  
  t<-c(rep(0, length(x15)),rep(2, length(x15)),rep(4, length(x15)),rep(6, length(x15)),rep(8, length(x15)),rep(10, length(x15)))
  amount<-c(x15,x30,x60,x120,x240,x360)
  tail<-rep(y,6)
  
  df <- data.frame(t,amount,tail)
  names(df) <- c('t', 'amount', 'tail')
  
  return(df)
}

#Definition of characteristics of the figure
w=1
size=24
wid_line=0.7
wid_clon=0.2
wid_comet=4
x11()


#Read experimental clonogenic data
survival_real_data<-read.table("../data/abrams_real_survival_3_6Gy.txt",header=F)
names(survival_real_data)<-c('dose','S','sd')


#Read experimental comet data
data<-matrix(data = scan("../data/abrams_comet_binning_40.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
t <- c(0, 2, 4, 6, 8, 10)
time<-c()
amount<-c()
mean_tail<-rep(seq(2,38, by=4),length(t))
for (i in 1:(length(t))) {
  time<-c(time,rep(t[i],length(data[i,])))
  amount<-c(amount,data[i,]*w)
}
comet_real_data<- data.frame(time,amount,mean_tail)
names(comet_real_data) <- c('t', 'amount', 'tail')


#Plot of clonogenic data from clonogenic calibration
survival_simu_data<-read.table("../analysis/synthetic/synthetic_survival_min_error_clon.txt", header=T)
survival_simu_data$fitting<-"Best~fit~according~to~epsilon[clonogenic]"
survival_simu_data$title<-"Clonogenic~assay~(survival)"
p1<-ggplot(survival_simu_data, aes(x=dose, y=S))+geom_line(size=wid_line, color="orange")+xlab("Dose / Gy")+ylab("Surviving fraction of cells")+scale_y_log10(breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1))+geom_errorbar(data=survival_real_data,aes(ymin=S-sd,ymax=S+sd),width=wid_clon,position=position_dodge(0.05),size=wid_line, color="blue")+geom_point(data=survival_real_data, color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(cols=vars(title),labeller = label_parsed)


#Plot of comet data from clonogenic calibration
comet_simu_data <- matrix(data = scan("../analysis/synthetic/synthetic_comet_min_error_clon.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
comet_simu_data <- comet_violin(comet_simu_data)
comet_simu_data$fitting<-"Best~fit~according~to~epsilon[clonogenic]"
comet_simu_data$title<-"Alkaline~comet~assay~(DNA~damage)"
p2<-ggplot(comet_simu_data, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line, color="orange")+xlab("Time / min")+ylab("Relative tail intensity / %")+geom_errorbar(data=comet_real_data,aes(xmin=t-amount,xmax=t+amount,y=tail),width=wid_comet,position=position_dodge(0.05),size=wid_line, color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(cols=vars(title),rows =vars(fitting),labeller = label_parsed)+scale_y_continuous(breaks=c(0,10,20,30,40))+scale_x_continuous(breaks=c(0,2,4,6,8,10),labels=c("15","30","60","120","240","360"))


#Plot of clonogenic data from comet calibration
survival_simu_data<-read.table("../analysis/synthetic/synthetic_survival_min_error_comet.txt", header=T)
survival_simu_data$fitting<-"Best~fit~according~to~epsilon[comet]"
survival_simu_data$title<-"Clonogenic~assay~(survival)"
p3<-ggplot(survival_simu_data, aes(x=dose, y=S))+geom_line(size=wid_line,color="orange")+xlab("Dose / Gy")+ylab("Surviving fraction of cells")+scale_y_log10(breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1))+geom_errorbar(data=survival_real_data,aes(ymin=S-sd,ymax=S+sd),width=wid_clon,position=position_dodge(0.05),size=wid_line,color="blue")+geom_point(data=survival_real_data,color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())


#Plot of comet data from comet calibration
comet_simu_data <- matrix(data = scan("../analysis/synthetic/synthetic_comet_min_error_comet.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
comet_simu_data <- comet_violin(comet_simu_data)
comet_simu_data$title<-"Alkaline~comet~assay~(DNA~damage)"
comet_simu_data$fitting<-"Best~fit~according~to~epsilon[comet]"
p4<-ggplot(comet_simu_data, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line,color="orange")+xlab("Time / min")+ylab("Relative tail intensity / %")+geom_errorbar(data=comet_real_data,aes(xmin=t-amount,xmax=t+amount,y=tail),width=wid_comet,position=position_dodge(0.05),size=wid_line,color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(rows =vars(fitting),labeller = label_parsed)+scale_y_continuous(breaks=c(0,10,20,30,40))+scale_x_continuous(breaks=c(0,2,4,6,8,10),labels=c("15","30","60","120","240","360"))


#Plot of clonogenic data from combined (eps=1) calibration
survival_simu_data<-read.table("../analysis/synthetic/synthetic_survival_min_error_combined_1.txt", header=T)
survival_simu_data$fitting<-"Best~fit~according~to~epsilon[combined]"
survival_simu_data$title<-"Clonogenic~assay~(survival)"
p5<-ggplot(survival_simu_data, aes(x=dose, y=S))+geom_line(size=wid_line,color="orange")+xlab("Dose / Gy")+ylab("Surviving fraction of cells")+scale_y_log10(breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1))+geom_errorbar(data=survival_real_data,aes(ymin=S-sd,ymax=S+sd),width=wid_clon,position=position_dodge(0.05),size=wid_line,color="blue")+geom_point(data=survival_real_data,color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())


#Plot of comet data from combined (eps=1) calibration
comet_simu_data <- matrix(data = scan("../analysis/synthetic/synthetic_comet_min_error_combined_1.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
comet_simu_data <- comet_violin(comet_simu_data)
comet_simu_data$title<-"Alkaline~comet~assay~(DNA~damage)"
comet_simu_data$fitting<-"Best~fit~according~epsilon[combined]~(xi==1)"
p6<-ggplot(comet_simu_data, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line,color="orange")+xlab("Time / min")+ylab("Relative tail intensity / %")+geom_errorbar(data=comet_real_data,aes(xmin=t-amount,xmax=t+amount,y=tail),width=wid_comet,position=position_dodge(0.05),size=wid_line,color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(rows =vars(fitting),labeller = label_parsed)+scale_y_continuous(breaks=c(0,10,20,30,40))+scale_x_continuous(breaks=c(0,2,4,6,8,10),labels=c("15","30","60","120","240","360"))


#Plot of clonogenic data from combined (eps=1/30) calibration
survival_simu_data<-read.table("../analysis/synthetic/synthetic_survival_min_error_combined_30.txt", header=T)
survival_simu_data$fitting<-"Best~fit~according~to~epsilon[combined]"
survival_simu_data$title<-"Clonogenic~assay~(survival)"
p7<-ggplot(survival_simu_data, aes(x=dose, y=S))+geom_line(size=wid_line,color="orange")+xlab("Dose / Gy")+ylab("Surviving fraction of cells")+scale_y_log10(breaks=c(0.4,0.5,0.6,0.7,0.8,0.9,1))+geom_errorbar(data=survival_real_data,aes(ymin=S-sd,ymax=S+sd),width=wid_clon,position=position_dodge(0.05),size=wid_line,color="blue")+geom_point(data=survival_real_data,color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())


#Plot of comet data from combined (eps=1/30) calibration
comet_simu_data <- matrix(data = scan("../analysis/synthetic/synthetic_comet_min_error_combined_30.txt", skip=0, sep=" ", fill = TRUE), ncol=10, byrow=TRUE)
comet_simu_data <- comet_violin(comet_simu_data)
comet_simu_data$title<-"Alkaline~comet~assay~(DNA~damage)"
comet_simu_data$fitting<-"Best~fit~according~to~epsilon[combined]~(xi==1/30)"
p8<-ggplot(comet_simu_data, aes(x=amount, y=tail, group=t))+geom_path(size=wid_line,color="orange")+xlab("Time / min")+ylab("Relative tail intensity / %")+geom_errorbar(data=comet_real_data,aes(xmin=t-amount,xmax=t+amount,y=tail),width=wid_comet,position=position_dodge(0.05),size=wid_line,color="blue")+theme_bw(base_size=size)+theme(strip.background = element_rect(colour="black", fill="white"),legend.position="top",legend.title=element_blank())+facet_grid(rows =vars(fitting),labeller = label_parsed)+scale_y_continuous(breaks=c(0,10,20,30,40))+scale_x_continuous(breaks=c(0,2,4,6,8,10),labels=c("15","30","60","120","240","360"))


#Save the figure
pdf("Fig_4.pdf",width=16,height=24)
ggarrange(p1+ xlab("")+theme(axis.text.x=element_blank()), p2 + xlab("")+theme(axis.text.x=element_blank()), p3+ xlab("")+theme(axis.text.x=element_blank()), p4+ xlab("")+theme(axis.text.x=element_blank()),p5+ xlab("")+theme(axis.text.x=element_blank()), p6+ xlab("")+theme(axis.text.x=element_blank()), p7, p8, labels = c("A", "B","C","D","E","F","G","H"), font.label = list(size = size),ncol = 2, nrow = 4)
graphics.off()

