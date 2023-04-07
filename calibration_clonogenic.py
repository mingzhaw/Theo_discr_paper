import sys
import numpy as np
import time
from mhr_model import mhr_model_clon
import logging
logging.basicConfig(level=logging.DEBUG)

# Define backend
#from abcpy.backends import BackendDummy as Backend
from abcpy.backends import BackendMPI as Backend
backend = Backend()

#Read the experimental clonogenic data
survival_data = np.loadtxt('data/abrams_real_survival_3_6Gy.txt',dtype=float)
data = [ np.log10(survival_data[:,1]) ]
datsize = len(data)

# Define Graphical Model
from abcpy.continuousmodels import Uniform, Normal
alpha = Uniform([[0.17], [ 2.00]], name='alpha')
cr    = Uniform([[0.00], [10.00]], name='cr')
ce    = Uniform([[0.00], [10.00]], name='ce')
mu_g  = Uniform([[0.00], [10.00]], name='mu_g')
gamma = Uniform([[0.00], [10.00]], name='gamma')
modelout = mhr_model_clon([datsize, alpha, cr, ce, mu_g, gamma], name = 'modelout')

# Define Statistics
from abcpy.statistics import Identity
statistics_calculator = Identity(degree=2, cross=False)

from abcpy.distances import Euclidean, LogReg
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
journal_sabc.save('sabc_journal_clon.jrnl')

print(" ")
print("Inference complete in %.3f seconds" % (time.time()-tsabc))
print("************************************")

#Save the calibration results
alpha_post_sample = np.array(journal_sabc.get_parameters()['alpha']).flatten()
cr_post_sample    = np.array(journal_sabc.get_parameters()['cr']).flatten()
ce_post_sample    = np.array(journal_sabc.get_parameters()['ce']).flatten()
mu_g_post_sample  = np.array(journal_sabc.get_parameters()['mu_g']).flatten()
gamma_post_sample = np.array(journal_sabc.get_parameters()['gamma']).flatten()

np.savetxt('results/sabc_calibration_clonogenic.dat',([alpha_post_sample, cr_post_sample, ce_post_sample, mu_g_post_sample, gamma_post_sample]))
