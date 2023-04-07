import sys
import numpy as np
import time
from mhr_model import mhr_model_combined_30
import logging
logging.basicConfig(level=logging.DEBUG)

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
data_clon = np.zeros(6)
data_clon[0] = np.log10(survival_data[0,1])
data_clon[1] = np.log10(survival_data[1,1]) 


#Merge the clonogenic and comet experimental data for a combined calibration
weight=1/np.sqrt(30)
data = np.zeros((6,11))
for i in range(0,6):
	for j in range(0,11):
		if j==10:
			data[i,j]=data_clon[i]
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


from abcpy.distances import Euclidean, LogReg,  Relative_mean
distance_calculator = Euclidean(statistics_calculator)


# Define kernel
from abcpy.perturbationkernel import DefaultKernel
kernel = DefaultKernel([alpha,cr,ce,mu_g,gamma])

# SABC ##
tsabc=time.time()
from abcpy.inferences import SABC
sampler = SABC([modelout], [distance_calculator], backend, kernel, seed = 1)
steps, epsilon, n_samples, n_samples_per_param = 250, np.array([0.5]), 5000, 1
journal_sabc = sampler.sample([data], steps, epsilon, n_samples, n_samples_per_param, ar_cutoff = 0.00001, full_output=0)
journal_sabc.save('sabc_journal_combined_30.jrnl')

print(" ")
print("Inference complete in %.3f seconds" % (time.time()-tsabc))
print("************************************")


#Save the calibration results
alpha_post_sample = np.array(journal_sabc.get_parameters()['alpha']).flatten()
cr_post_sample    = np.array(journal_sabc.get_parameters()['cr']).flatten()
ce_post_sample    = np.array(journal_sabc.get_parameters()['ce']).flatten()
mu_g_post_sample  = np.array(journal_sabc.get_parameters()['mu_g']).flatten()
gamma_post_sample = np.array(journal_sabc.get_parameters()['gamma']).flatten()

np.savetxt('results/sabc_calibration_combined_30.dat',([alpha_post_sample, cr_post_sample, ce_post_sample, mu_g_post_sample, gamma_post_sample]))

