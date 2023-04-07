#Open required libraries
library(reshape2)
library(ggplot2)
library(ggh4x)


#Make sure that there is no open graphs
graphics.off()


#Function to read the data
load_file <- function(file, sigma)
{
  raw_data<-read.table(file)
  raw_data<-t(unname(raw_data))
  raw_data<-data.frame(raw_data)
  names(raw_data)<-c('alpha~(Gy^-1)','c[r]~(h^-1)','c[e]~(h^-1)','mu[Gamma]~(Gy^-1)','gamma~(h^-1)')
  raw_data<-melt(raw_data[,3])
  raw_data$fitting<-sigma
  return(raw_data)
}


#Read the results obtained from the simulations
data_initial <- load_file("../additional_simulations/sabc_calibration_combined_30_initial.dat","Initial~distribution")
data_final <- load_file("../results/sabc_calibration_combined_30.dat","Final~distribution")
data<-rbind(data_initial,data_final)
data$fitting <- factor(data$fitting,levels=c("Initial~distribution","Final~distribution"))


#Produce and save the figure
pdf("Fig_3B.pdf",width=6,height=12)
print(
  ggplot(data, aes(x=value, group=fitting))+ 
    geom_histogram(aes(y = after_stat(..density..*0.1), group=fitting),breaks = seq(0,10,0.1))+
    theme_bw(base_size=23)+
    xlab(expression(c[e]~(h^-1)))+
    ylab("Relative Frequency")+
    ylim(c(0,0.45))+
    facet_grid(rows = vars(fitting),scales="free", labeller = label_parsed)+
    scale_x_continuous(breaks=c(0,2,4,6,8,10),labels=c("0","2","4","6","8","10"))
    
)
graphics.off()
