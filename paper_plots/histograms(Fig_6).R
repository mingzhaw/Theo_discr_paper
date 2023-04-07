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
  raw_data<-melt(raw_data[,1:5])
  raw_data$fitting<-sigma
  return(raw_data)
}


#Read the results obtained from the simulations
data_weyland <- load_file("weyland_combined_parameters.dat","Weyland~et~al.~fitting")
data_combined <- load_file("../results/sabc_calibration_combined_30.dat","Combined~fitting~(xi==1/30)")
data_discriminators <- load_file("../results/sabc_after_dr_discriminators_final.dat","Theoretical~discriminators")
data<-rbind(data_weyland,data_combined,data_discriminators)
data$fitting <- factor(data$fitting,levels=c("Weyland~et~al.~fitting","Combined~fitting~(xi==1/30)","Theoretical~discriminators"))


#Produce and save the figure
pdf("Fig_6.pdf",width=20,height=12)
print(
  ggplot(data, aes(x=value, group=fitting))+ 
    geom_histogram(data = data[data$variable == "alpha~(Gy^-1)",], aes(y = after_stat(..density..*0.05), group=fitting),breaks = seq(0,2,0.05)) +
    geom_histogram(data = data[data$variable == "c[r]~(h^-1)",], aes(y = after_stat(..density..*0.5), group=fitting),breaks = seq(0,10,0.5))+
    geom_histogram(data = data[data$variable == "c[e]~(h^-1)",], aes(y = after_stat(..density..*0.1), group=fitting),breaks = seq(0,10,0.1))+
    geom_histogram(data = data[data$variable == "mu[Gamma]~(Gy^-1)",], aes(y = after_stat(..density..*0.2), group=fitting),breaks = seq(0,10,0.2)) +
    geom_histogram(data = data[data$variable == "gamma~(h^-1)",], aes(y = after_stat(..density..*0.5), group=fitting),breaks = seq(0,10,0.5))+
    theme_bw(base_size=23)+
    xlab("Parameter range")+
    ylab("Relative Frequency")+
    facet_grid(cols=vars(variable),rows = vars(fitting),scales="free", labeller = label_parsed)+
    facetted_pos_scales(x = list(variable != "alpha~(Gy^-1)" ~ scale_x_continuous(breaks=c(0,2,4,6,8,10),labels=c("0","2","4","6","8","10"))))
)
graphics.off()