from scipy.integrate import odeint
import numpy as np
from math import exp, log10, log


#Functions to get the synthetic comets 
def get_index_from_minutes(minutes):
    return np.abs(ts-minutes/60.).argmin()


def get_normalized_row(Ps, index):
    row = Ps[index,0:10]
    return row/np.sum(row)


#Definition of the MHR model
stop_time = 10  # in h (We need 6h for comet but >= 10 hours for survival limits)
dose_rate = 360 # in Gy/h


def R(time, dose, dose_rate):
    rad_time = dose / dose_rate
    if time >= 0 and time < rad_time:
        return dose/rad_time
    else:
        return 0.

def dP_dt(P, t, alpha, cr, ce, mu_g, gamma, dose, dose_rate):
    (N, L01, L02, L03, L04, L05, L06, L07, L08, L09, Gamma) = P
    r_GL = cr * exp(-mu_g*Gamma)
    R_t = R(t, dose, dose_rate)
    
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
P0 = [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]




###################################################
#                                                 #
#   Synthetic comets for different alpha values   #
#                                                 #
###################################################


#Synthetic comet produced after irradiation for alpha = 0.3 Gy-1
fun=lambda t, x: dP_dt(t,x,0.3,0.,0.,0.,0.,6.,360)
Ps = odeint(fun, P0, ts, hmax = 1./3600 )
comet03=get_normalized_row(Ps, get_index_from_minutes(15))


#Synthetic comet produced after irradiation for alpha = 1.0 Gy-1
fun=lambda t, x: dP_dt(t,x,1.,0.,0.,0.,0.,6.,360)
Ps = odeint(fun, P0, ts, hmax = 1./3600 )
comet10=get_normalized_row(Ps, get_index_from_minutes(15))


#Synthetic comet produced after irradiation for alpha = 1.5 Gy-1
fun=lambda t, x: dP_dt(t,x,1.5,0.,0.,0.,0.,6.,360)
Ps = odeint(fun, P0, ts, hmax = 1./3600 )
comet15=get_normalized_row(Ps, get_index_from_minutes(15))


#Synthetic comet produced after irradiation for alpha = 2.0 Gy-1
fun=lambda t, x: dP_dt(t,x,2.,0.,0.,0.,0.,6.,360)
Ps = odeint(fun, P0, ts, hmax = 1./3600 )
comet20=get_normalized_row(Ps, get_index_from_minutes(15))


#Save results
np.savetxt('comet_alpha.dat',([comet03, comet10, comet15, comet20]))




#############################################################
#                                                           #
#   Generation of survival curves at different dose-rates   #
#                                                           #
#############################################################


#Read a parameter set which passes the theoretical discriminators as example
data= np.loadtxt('../results/sabc_after_dr_discriminators_final.dat', dtype='float')

alpha = data[0,10]
cr = data[1,10]
ce = data[2,10]
mu_g = data[3,10]
gamma = data[4,10]


#Definition of vectors needed for simulating and saving the survival curves
D =np.linspace(0.5,8, 16)
D1 =np.linspace(0.,8, 17)

survival_001 =np.zeros(len(D)+1)
survival_01 =np.zeros(len(D)+1)
survival_2 =np.zeros(len(D)+1)
survival_20 =np.zeros(len(D)+1)


survival_001[0]=1
survival_01[0]=1
survival_2[0]=1
survival_20[0]=1


#Loop for generating the survival curves
for i in range(0,len(D)):
    fun1=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,D[i],0.1*60)
    Ps1 = odeint(fun1, P0, ts, hmax = 1./3600 )
    fun2=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,D[i],2*60)
    Ps2 = odeint(fun2, P0, ts, hmax = 1./3600 )
    fun3=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,D[i],20*60)
    Ps3 = odeint(fun3, P0, ts, hmax = 1./3600 )
       
    survival_01[i+1]=Ps1[-1,0]
    survival_2[i+1]=Ps2[-1,0]
    survival_20[i+1]=Ps3[-1,0]
       
    dose_rate = 0.01*60 # in Gy/h
    rad_time=D[i]/dose_rate
    ts1 = np.linspace(0, stop_time+rad_time, round(stop_time+rad_time)*60*60)
    fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,D[i], dose_rate)
    Ps = odeint(fun, P0, ts1, hmax = 1./3600 )
    survival_001[i+1] = Ps[-1,0]


#Save the results
np.savetxt('survival_after_dr_discriminators.txt',np.transpose([D1, survival_01, survival_2, survival_20, survival_001]), header="dose S_01 S_2 S_20 S_001", comments="")
