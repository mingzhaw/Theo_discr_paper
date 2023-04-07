# Running this program requires the selection of mymodel_comet
# in mhr_model.py!

import sys
import numpy as np
import time
from mhr_model import mhr_model_comet
import logging
logging.basicConfig(level=logging.DEBUG)

# Define backend
#from abcpy.backends import BackendDummy as Backend
from abcpy.backends import BackendMPI as Backend
backend = Backend()

#Read the experimental comet data
max_damage_mapping = int(40)

data = np.zeros((6,10))
i = 0
for t in ('15','30','60','120','240','360'):
    row = np.loadtxt('data/abrams_real_2_%s.txt' % t)
    hist = np.histogram(row, bins=10, range=(0,max_damage_mapping))[0]
    data[i,] = hist/np.sum(hist)
    i += 1
data = [ data ]
datsize = len(data)

# Define Graphical Model
from abcpy.continuousmodels import Uniform, Normal
alpha = Uniform([[0.17], [ 2.00]], name='alpha')
cr    = Uniform([[0.00], [10.00]], name='cr')
ce    = Uniform([[0.00], [10.00]], name='ce')
mu_g  = Uniform([[0.00], [10.00]], name='mu_g')
gamma = Uniform([[0.00], [10.00]], name='gamma')
modelout = mhr_model_comet([datsize, alpha, cr, ce, mu_g, gamma], name = 'modelout')

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
journal_sabc.save('sabc_journal_comet.jrnl')

print(" ")
print("Inference complete in %.3f seconds" % (time.time()-tsabc))
print("************************************")


#Save the calibration results
alpha_post_sample = np.array(journal_sabc.get_parameters()['alpha']).flatten()
cr_post_sample    = np.array(journal_sabc.get_parameters()['cr']).flatten()
ce_post_sample    = np.array(journal_sabc.get_parameters()['ce']).flatten()
mu_g_post_sample  = np.array(journal_sabc.get_parameters()['mu_g']).flatten()
gamma_post_sample = np.array(journal_sabc.get_parameters()['gamma']).flatten()

np.savetxt('results/sabc_calibration_comet.dat', ([alpha_post_sample, cr_post_sample, ce_post_sample, mu_g_post_sample, gamma_post_sample]))

