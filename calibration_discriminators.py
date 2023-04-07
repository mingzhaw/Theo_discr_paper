import sys
import numpy as np
import time
from mhr_model import mhr_model_combined_30
import logging
logging.basicConfig(level=logging.DEBUG)
from scipy.integrate import odeint
from math import exp, log10, log
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# Define backend
#from abcpy.backends import BackendDummy as Backend
from abcpy.backends import BackendMPI as Backend
backend = Backend()


#Read the experimental comet data
max_damage_mapping = int(40)
data_comet = np.zeros((6,10))
i = 0
for t in ('15','30','60','120','240','360'):
    row = np.loadtxt('data/abrams_real_2_%s.txt' % t)
    hist = np.histogram(row, bins=10, range=(0,max_damage_mapping))[0]
    data_comet[i,] = hist/np.sum(hist)
    i += 1


#Read the experimental clonogenic data
survival_data = np.loadtxt('data/abrams_real_survival_3_6Gy.txt',dtype=float)
data_clon_log = np.zeros(6)
data_clon_log[0] = np.log10(survival_data[0,1])
data_clon_log[1] = np.log10(survival_data[1,1]) 
data_clon =  survival_data[:,1] 
error_clon =  survival_data[:,2] 


#Merge the clonogenic and comet experimental data for a combined calibration
weight=1/np.sqrt(30)
data = np.zeros((6,11))
for i in range(0,6):
	for j in range(0,11):
		if j==10:
			data[i,j]=data_clon_log[i]
		else:
			data[i,j]=weight*data_comet[i,j]

data = [ data ]
datsize = len(data)


# Define Graphical Model
from abcpy.continuousmodels import Uniform, Normal
alpha = Uniform([[0.17], [ 2.00]], name='alpha')
cr    = Uniform([[0.00], [10.00]], name='cr')
ce    = Uniform([[0.00], [10.00]], name='ce')
mu_g  = Uniform([[0.00], [10.00]], name='mu_g')
gamma = Uniform([[0.00], [10.00]], name='gamma')
modelout = mhr_model_combined_30([datsize, alpha, cr, ce, mu_g, gamma], name = 'modelout')

# Define Statistics
from abcpy.statistics import Identity
statistics_calculator = Identity(degree=2, cross=False)


from abcpy.distances import Euclidean, LogReg
#distance_calculator = LogReg(statistics_calculator)
distance_calculator = Euclidean(statistics_calculator)


# Define kernel
from abcpy.perturbationkernel import DefaultKernel
kernel = DefaultKernel([alpha,cr,ce,mu_g,gamma])

# SABC ##
tsabc=time.time()
from abcpy.inferences import SABC



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



#Linear Quadratic (LQ) model for the Additional discriminator
def LQ(x, a, b):

    return -a*x -b*x*x
    
#Linear model for the low dose-rate discriminator
model = LinearRegression(fit_intercept=False)


#Definition of the calibration parameters    
steps, epsilon, n_samples, n_samples_per_param = 250, np.array([0.5]), 200, 1

#Definition of counter (j) of accepted parameters sets and number of parameters (N) at which the calibration stops
j=0
N=50


#Initialization of vectors which pass the different theoretical discriminators
alpha_post_sample = np.zeros(N)
cr_post_sample    = np.zeros(N)
ce_post_sample    = np.zeros(N)
mu_g_post_sample  = np.zeros(N)
gamma_post_sample = np.zeros(N)



while j<N:
    sampler = SABC([modelout], [distance_calculator], backend, kernel)
    journal_sabc = sampler.sample([data], steps, epsilon, n_samples, n_samples_per_param, ar_cutoff = 0.00001, full_output=0)
    
    alpha = np.array(journal_sabc.get_parameters()['alpha']).flatten()
    cr    = np.array(journal_sabc.get_parameters()['cr']).flatten()
    ce    = np.array(journal_sabc.get_parameters()['ce']).flatten()
    mu_g  = np.array(journal_sabc.get_parameters()['mu_g']).flatten()
    gamma = np.array(journal_sabc.get_parameters()['gamma']).flatten()
    
    print("accepted sets=",j)
    
    for i in range(0,n_samples):
        
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
        
        #Conditions I-IV    
        if survival_2[0]>=survival_20[0] and survival_2[1]>=survival_20[1] and survival_2[2]>=survival_20[2] and survival_01[0]>survival_2[0] and survival_01[1]>survival_2[1] and survival_01[2]>survival_2[2] and (survival_01[2]-survival_2[2])>(survival_01[1]-survival_2[1]) and (survival_01[1]-survival_2[1])>(survival_01[0]-survival_2[0]):
                        
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
            
                #Additional discriminator
                for dose in range(1,11):
                    fun=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],dose,360, Fract=False)
                    Ps = odeint(fun, P0, ts, hmax = 1./3600 )
                    survival_log[0,dose] = np.log(Ps[-1,0])
                    survival[dose] = Ps[-1,0]
        
                y = [survival_log[0,0], survival_log[0,1], survival_log[0,2], survival_log[0,3], survival_log[0,4], survival_log[0,5], survival_log[0,6], survival_log[0,7], survival_log[0,8], survival_log[0,9], survival_log[0,10]]
        
                popt, pcov = curve_fit(LQ, D_LQ, y)
                y_pred = LQ(D_LQ, *popt)
    
                #Conditions VII and VIII    
                if abs(data_clon[0]-survival[3])<=error_clon[0] and abs(data_clon[1]-survival[6])<=error_clon[1] and r2_score(y, y_pred)>=0.99 and  popt[0]>0 and popt[1]>0:
    
    
                    #Fractionation discriminator    
                    fun=lambda t, x: dP_dt(t,x,alpha[i],cr[i],ce[i],mu_g[i],gamma[i],6., 360., Fract=True)
                    Ps = odeint(fun, P0, ts1, hmax = 1./3600 )
                    survival_fract = Ps[-1,0]
                
                    #Condition VI
                    if survival_fract>survival[6]:
                        if j<N:
                            alpha_post_sample[j] = alpha[i]
                            cr_post_sample[j]    = cr[i]
                            ce_post_sample[j]    = ce[i]
                            mu_g_post_sample[j]  = mu_g[i]
                            gamma_post_sample[j] = gamma[i]
                            j+=1
                        else:
                            alpha_post_sample = np.append(alpha_post_sample,alpha[i])
                            cr_post_sample    = np.append(cr_post_sample,cr[i])
                            ce_post_sample    = np.append(ce_post_sample,ce[i])
                            mu_g_post_sample  = np.append(mu_g_post_sample,mu_g[i])
                            gamma_post_sample = np.append(gamma_post_sample,gamma[i])
                            j+=1
    
    
    print("accepted sets=",j)
    del(journal_sabc)
    print(" ")
    print("Inference complete in %.3f seconds" % (time.time()-tsabc))
    print("************************************")
    
    np.savetxt('results/sabc_calibration_discriminators.dat',([alpha_post_sample, cr_post_sample, ce_post_sample, mu_g_post_sample, gamma_post_sample]))
    
    
    

