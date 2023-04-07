import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from scipy.integrate import odeint
from math import exp, log10, log


#Functions to get the synthetic comets 
def get_index_from_minutes(minutes):
    return np.abs(ts-minutes/60.).argmin()


def get_normalized_row(Ps, index):
    row = Ps[index,0:10]
    return row/np.sum(row)


#Linear Quadratic (LQ) model for the Additional discriminator
def LQ(x, a, b):

    return -a*x -b*x*x


#Linear model for the low dose-rate discriminator
model = LinearRegression(fit_intercept=False)


#Read the experimental comet data
max_damage_mapping = int(40)
data_comet = np.zeros((6,10))
i = 0
for t in ('15','30','60','120','240','360'):
    row = np.loadtxt('../data/abrams_real_2_%s.txt' % t)
    hist = np.histogram(row, bins=10, range=(0,max_damage_mapping))[0]
    data_comet[i,] = hist/np.sum(hist)
    i += 1


#Read the experimental clonogenic data
survival_data = np.loadtxt('../data/abrams_real_survival_3_6Gy.txt',dtype=float)
data_clon =  survival_data[:,1] 
data_clon_error =  survival_data[:,2] 



#Definition of the MHR model
stop_time = 10  # in h (We need 6h for comet but >= 10 hours for survival limits)
dose_rate = 360 # in Gy/h


def R1(time, dose, dose_rate):
    rad_time = dose / dose_rate
    if time >= 0 and time < rad_time:
        return dose/rad_time
    else:
        return 0.

def R2(time, dose, dose_rate):
    rad_time = dose / dose_rate
    if time >= 0 and time < rad_time/2:
        return dose/rad_time
    elif time >= (rad_time/2 + 2) and time < (rad_time + 2):
    	return dose/rad_time
    else:
        return 0.

def dP_dt(P, t, alpha, cr, ce, mu_g, gamma, dose, dose_rate, Fract=False):
    (N, L01, L02, L03, L04, L05, L06, L07, L08, L09, Gamma) = P
    r_GL = cr * exp(-mu_g*Gamma)
    
    if(Fract):
        R_t = R2(t, dose, dose_rate)
    else:
        R_t = R1(t, dose, dose_rate)
    
    return [
                           - alpha*R_t*N                       + r_GL*L01,
             alpha*R_t*N   - alpha*R_t*L01 - r_GL*L01 - ce*L01 + r_GL*L02,
             alpha*R_t*L01 - alpha*R_t*L02 - r_GL*L02 - ce*L02 + r_GL*L03,
             alpha*R_t*L02 - alpha*R_t*L03 - r_GL*L03 - ce*L03 + r_GL*L04,
             alpha*R_t*L03 - alpha*R_t*L04 - r_GL*L04 - ce*L04 + r_GL*L05,
             alpha*R_t*L04 - alpha*R_t*L05 - r_GL*L05 - ce*L05 + r_GL*L06,
             alpha*R_t*L05 - alpha*R_t*L06 - r_GL*L06 - ce*L06 + r_GL*L07,
             alpha*R_t*L06 - alpha*R_t*L07 - r_GL*L07 - ce*L07 + r_GL*L08,
             alpha*R_t*L07 - alpha*R_t*L08 - r_GL*L08 - ce*L08 + r_GL*L09,
             alpha*R_t*L08 - alpha*R_t*L09 - r_GL*L09 - ce*L09,
             R_t - gamma * Gamma
           ]

ts = np.linspace(0, stop_time, stop_time*60*60)
ts1 = np.linspace(0, stop_time+2, (stop_time+2)*60*60)
P0 = [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]


#Read the parameter sets obtained from the model calibration
data= np.loadtxt('../results/sabc_calibration_comet.dat', dtype='float')

alpha = data[0,:]
cr = data[1,:]
ce = data[2,:]
mu_g = data[3,:]
gamma = data[4,:]

N=len(alpha)



#Definition of vectors / variables needed for applying the discriminators
survival_log = np.zeros((1,11))
survival_log[0,0] = np.log(1)
D = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]).reshape((-1, 1))
survival = np.zeros(11)
survival_01 =np.zeros(3)
survival_2 =np.zeros(3)
survival_20 =np.zeros(3)
D1=[3,6,9]
D_LQ = np.linspace(0, 10, 11)



#Initialization of the counters to track if the parameter sets pass the different theoretical discriminators
survival_dr_discr   = np.zeros(N)
low_dr_discr        = np.zeros(N)
fractionation_discr = np.zeros(N)
additional_discr    = np.zeros(N)
all_discr   = np.zeros(N)


#Initialization of vectors to store the errors 
error_clon  = np.zeros(N)
error_comet = np.zeros(N)
    

#Loop for passing all the parameter sets through the theoretical discriminators
for i in range(0,N):
    print("Parameter set %i"%(i+1)+" out of %i" %(N), end="\r") 
    
    #Survival at different dose-rates discriminator
    k=0
    for dose in D1:
        fun1=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],dose,0.1*60, Fract=False)
        Ps1 = odeint(fun1, P0, ts, hmax = 1./3600 )
        fun2=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],dose,2*60, Fract=False)
        Ps2 = odeint(fun2, P0, ts, hmax = 1./3600 )
        fun3=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],dose,20*60, Fract=False)
        Ps3 = odeint(fun3, P0, ts, hmax = 1./3600 )
        
        survival_01[k]=Ps1[-1,0]
        survival_2[k]=Ps2[-1,0]
        survival_20[k]=Ps3[-1,0]
        k+=1
        
    #Condition I    
    if survival_2[0]>=survival_20[0] and survival_2[1]>=survival_20[1] and survival_2[2]>=survival_20[2]:
        #Condition II
        if survival_01[0]>survival_2[0] and survival_01[1]>survival_2[1] and survival_01[2]>survival_2[2]:
            #Conditions III and IV
            if (survival_01[2]-survival_2[2])>(survival_01[1]-survival_2[1]) and (survival_01[1]-survival_2[1])>(survival_01[0]-survival_2[0]):
                survival_dr_discr[i] = 1                
        
        
        
    #Low dose-rate discriminator
    for dose in range(1,11):
        rad_time=dose/(0.01*60)
        ts3 = np.linspace(0, stop_time+rad_time, round(stop_time+rad_time)*60*60)
        fun=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],dose,0.01*60, Fract=False)
        Ps = odeint(fun, P0, ts3, hmax = 1./3600 )
        survival_log[0,dose] = np.log(Ps[-1,0])
    
    y = np.array([survival_log[0,0], survival_log[0,1], survival_log[0,2], survival_log[0,3], survival_log[0,4], survival_log[0,5], survival_log[0,6], survival_log[0,7], survival_log[0,8], survival_log[0,9], survival_log[0,10]])
    model.fit(D, y)
    
    #Condition V    
    if model.score(D, y)>=0.99:
        low_dr_discr[i] = 1
        
        
            
    #Additional discriminator
    for dose in range(1,11):
        fun=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],dose,360, Fract=False)
        Ps = odeint(fun, P0, ts, hmax = 1./3600 )
        survival_log[0,dose] = np.log(Ps[-1,0])
        survival[dose] = Ps[-1,0]
        
        #Produce synthetic comets
        if dose==6:
            comet15=get_normalized_row(Ps, get_index_from_minutes( 15))
            comet30=get_normalized_row(Ps, get_index_from_minutes( 30))
            comet60=get_normalized_row(Ps, get_index_from_minutes( 60))
            comet120=get_normalized_row(Ps, get_index_from_minutes(120))
            comet240=get_normalized_row(Ps, get_index_from_minutes(240))
            comet360=get_normalized_row(Ps, get_index_from_minutes(360))
                
    y = [survival_log[0,0], survival_log[0,1], survival_log[0,2], survival_log[0,3], survival_log[0,4], survival_log[0,5], survival_log[0,6], survival_log[0,7], survival_log[0,8], survival_log[0,9], survival_log[0,10]]
        
    popt, pcov = curve_fit(LQ, D_LQ, y)
    y_pred = LQ(D_LQ, *popt)
    
    #Condition VII    
    if abs(data_clon[0]-survival[3])<=data_clon_error[0] and abs(data_clon[1]-survival[6])<=data_clon_error[1]: 
        #Condition VIII
        if r2_score(y, y_pred)>=0.99 and popt[0]>0 and popt[1]>0:
            additional_discr[i] = 1  
    
  
    
    
    #Fractionation discriminator    
    fun=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],6., 360., Fract=True)
    Ps = odeint(fun, P0, ts1, hmax = 1./3600 )
    survival_fract = Ps[-1,0]
    
    #Condition VI
    if survival_fract>survival[6]:
        fractionation_discr[i] = 1
    
    
    #Check if the parameter set passes all the discriminators
    if survival_dr_discr[i]==1 and low_dr_discr[i]==1 and fractionation_discr[i]==1 and additional_discr[i]==1:
        all_discr[i] = 1
                       
    
    
    #Calculate comet and clonogenic errors        
    for j in range(0,10):
        error_comet[i]=error_comet[i]+(comet15[j]-data_comet[0,j])**2+(comet30[j]-data_comet[1,j])**2+(comet60[j]-data_comet[2,j])**2+(comet120[j]-data_comet[3,j])**2+(comet240[j]-data_comet[4,j])**2+(comet360[j]-data_comet[5,j])**2
    
    
    error_clon[i]=(np.log10(data_clon[0])-np.log10(survival[3]))**2+(np.log10(data_clon[1])-np.log10(survival[6]))**2    



#Save the parameters sets which pases the all the discriminators
if np.sum(all_discr)>0:
    np.savetxt('after_discriminators_comet.dat', ([alpha[all_discr==1], cr[all_discr==1], ce[all_discr==1], mu_g[all_discr==1], gamma[all_discr==1]]))


#Save results of the analysis
file = open("analysis_comet.txt", "w")

file.write("Percentage of accepted paremeters after Survival at different dose-rates discriminator: %.3f" %(np.sum(survival_dr_discr)*100/N)  + "\n")
file.write("Percentage of accepted paremeters after Low dose-rate discriminator: %.3f"  %(np.sum(low_dr_discr)*100/N)  + "\n")
file.write("Percentage of accepted paremeters after Fractionation discriminator: %.3f"  %(np.sum(fractionation_discr)*100/N)  + "\n")
file.write("Percentage of accepted paremeters after Additional discriminator: %.3f"  %(np.sum(additional_discr)*100/N)  + "\n")
file.write("Percentage of accepted paremeters after all discriminators: %.3f" %(np.sum(all_discr)*100/N)  + "\n")
file.write("\n")
file.write("Clonogenic error (minimum, maximum, mean): %.3e" %(np.min(error_clon)) + "  %.3e"   %(np.max(error_clon)) + "  %.3e"  %(np.mean(error_clon)) + "  "  + "\n")
file.write("Comet error (minimum, maximum, mean): %.3e" %(np.min(error_comet)) + "  %.3e"   %(np.max(error_comet)) + "  %.3e"  %(np.mean(error_comet)) + "  "  + "\n")
    
file.close()                       



#Save information of the parameter set which provide a smaller error
np.savetxt('parameters/parameters_min_error_comet.dat',([alpha[np.argmin(error_comet)], cr[np.argmin(error_comet)], ce[np.argmin(error_comet)], mu_g[np.argmin(error_comet)], gamma[np.argmin(error_comet)], error_clon[np.argmin(error_comet)], error_comet[np.argmin(error_comet)]]))   

alpha = alpha[np.argmin(error_comet)]
cr    = cr[np.argmin(error_comet)]
ce    = ce[np.argmin(error_comet)]
mu_g  = mu_g[np.argmin(error_comet)]
gamma = gamma[np.argmin(error_comet)]


D =np.linspace(0.,8, 17)
survival = np.zeros(17)

    
for j in range(0,len(D)):
    fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,D[j],360, Fract=False)
    Ps = odeint(fun, P0, ts, hmax = 1./3600 )
    
    survival[j] = Ps[-1,0]
        
    if D[j]==6.:
        comet15=get_normalized_row(Ps, get_index_from_minutes( 15))
        comet30=get_normalized_row(Ps, get_index_from_minutes( 30))
        comet60=get_normalized_row(Ps, get_index_from_minutes( 60))
        comet120=get_normalized_row(Ps, get_index_from_minutes(120))
        comet240=get_normalized_row(Ps, get_index_from_minutes(240))
        comet360=get_normalized_row(Ps, get_index_from_minutes(360))

        np.savetxt('synthetic/synthetic_comet_min_error_comet.txt',([comet15, comet30, comet60, comet120, comet240, comet360]))
                       
np.savetxt('synthetic/synthetic_survival_min_error_comet.txt',np.transpose([D, survival]), header="dose S", comments="")                             
