from abcpy.probabilisticmodels import ProbabilisticModel, Continuous, InputConnector
import numpy as np
from scipy.integrate import odeint
from math import exp, log10, log

stop_time = 10  # in h (We need 6h for comet but >= 10 hours for survival limits)
dose_rate = 360 # in Gy/h
time_gap  = 10/60  # Duration of time gap, in hours


def R(time, dose):
    rad_time = dose / dose_rate
    if time >= 0 and time < rad_time:
        return dose/rad_time
    else:
        return 0.
        

def dP_dt(P, t, alpha, cr, ce, mu_g, gamma, dose):
    (N, L01, L02, L03, L04, L05, L06, L07, L08, L09, Gamma) = P
    R_t = R(t, dose)
    r_GL = cr * exp(-mu_g*Gamma)
    
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
            

def get_index_from_minutes(minutes):
    return np.abs(ts-minutes/60.).argmin()


def get_normalized_row(Ps, index):
    row = Ps[index,0:10]
    return row/np.sum(row)




class mhr_model_clon(ProbabilisticModel, Continuous):
    """
    This class is an re-implementation of the `abcpy.continousmodels.Normal` for documentation purposes.
    """

    def __init__(self, parameters, name='IterMap'):
        # We expect input of type parameters = [alpha,cr,ce,mu_g,gamma]
        if not isinstance(parameters, list):
            raise TypeError('Input of myModel is of type list')

        if len(parameters) != 6:
            raise RuntimeError('Input list must be of length 6, containing [datsize, alpha, cr, ce, mu_g, gamma].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 6:
            raise ValueError('Number of parameters of myModel must be 6.')

        # Check whether input is from correct domain (FIXME)
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]

        if not isinstance(datsize, np.int64):
            return False
        if datsize <= 0 or alpha < 0 or cr < 0 or ce < 0 or mu_g < 0 or gamma < 0:
            return False

        return True

    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]
        
        # Do the actual forward simulation
        vector_of_k_samples = self.mymodel_clonogenic(datsize, alpha, cr, ce, mu_g, gamma, k) # For clonogenic simulation
	# Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def mymodel_clonogenic(self, datsize, alpha, cr, ce, mu_g, gamma, k):

        # import warnings
        # warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

        result = []
        for ind_k in range(k):
            survival = []
            
            for dose in (3,6):
                fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,dose)
                Ps = odeint(fun, P0, ts, hmax = 1./3600 )
                survival.append(np.log10(Ps[-1,0]))

            result.append( survival )
        return result


  
class mhr_model_comet(ProbabilisticModel, Continuous):
    """
    This class is an re-implementation of the `abcpy.continousmodels.Normal` for documentation purposes.
    """

    def __init__(self, parameters, name='IterMap'):
        # We expect input of type parameters = [alpha,cr,ce,mu_g,gamma]
        if not isinstance(parameters, list):
            raise TypeError('Input of myModel is of type list')

        if len(parameters) != 6:
            raise RuntimeError('Input list must be of length 6, containing [datsize, alpha, cr, ce, mu_g, gamma].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 6:
            raise ValueError('Number of parameters of myModel must be 6.')

        # Check whether input is from correct domain (FIXME)
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]

        if not isinstance(datsize, np.int64):
            return False
        if datsize <= 0 or alpha < 0 or cr < 0 or ce < 0 or mu_g < 0 or gamma < 0:
            return False

        return True

    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]
        
        # Do the actual forward simulation
        vector_of_k_samples = self.mymodel_comet(datsize, alpha, cr, ce, mu_g, gamma, k) # For comet simulation
        # Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    def mymodel_comet(self, datsize, alpha, cr, ce, mu_g, gamma, k):

        # import warnings
        # warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

        result = []
        for ind_k in range(k):
            # Simulate at 6 Gy
            fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,6)
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            result.append( (
                get_normalized_row(Ps, get_index_from_minutes( 15) ),
                get_normalized_row(Ps, get_index_from_minutes( 30) ),
                get_normalized_row(Ps, get_index_from_minutes( 60) ),
                get_normalized_row(Ps, get_index_from_minutes(120) ),
                get_normalized_row(Ps, get_index_from_minutes(240) ),
                get_normalized_row(Ps, get_index_from_minutes(360) )
                )
            )
        return result

        
class mhr_model_combined_1(ProbabilisticModel, Continuous):
    """
    This class is an re-implementation of the `abcpy.continousmodels.Normal` for documentation purposes.
    """

    def __init__(self, parameters, name='IterMap'):
        # We expect input of type parameters = [alpha,cr,ce,mu_g,gamma]
        if not isinstance(parameters, list):
            raise TypeError('Input of myModel is of type list')

        if len(parameters) != 6:
            raise RuntimeError('Input list must be of length 6, containing [datsize, alpha, cr, ce, mu_g, gamma].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 6:
            raise ValueError('Number of parameters of myModel must be 6.')

        # Check whether input is from correct domain (FIXME)
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]

        if not isinstance(datsize, np.int64):
            return False
        if datsize <= 0 or alpha < 0 or cr < 0 or ce < 0 or mu_g < 0 or gamma < 0:
            return False

        return True

    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]
        
        # Do the actual forward simulation
        vector_of_k_samples = self.mymodel_combined(datsize, alpha, cr, ce, mu_g, gamma, k) # For combined simulation
	# Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    

    def mymodel_combined(self, datsize, alpha, cr, ce, mu_g, gamma, k):

        # import warnings
        # warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
        
        weight = 1
        result = np.zeros((6,11))
        for ind_k in range(k):
            #Simulate at 3 Gy (clonogenic)
            fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,3.)
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            result[0,10] = np.log10(Ps[-1,0])
            
            # Simulate at 6 Gy (comet and clonogenic and comet)
            fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,6.)
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            result[1,10] = np.log10(Ps[-1,0])
            
            result[0,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes( 15) )
            result[1,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes( 30) )
            result[2,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes( 60) )
            result[3,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes(120) )
            result[4,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes(240) )
            result[5,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes(360) )
                
            result = [result]
	
        return result
        
        

class mhr_model_combined_30(ProbabilisticModel, Continuous):
    """
    This class is an re-implementation of the `abcpy.continousmodels.Normal` for documentation purposes.
    """

    def __init__(self, parameters, name='IterMap'):
        # We expect input of type parameters = [alpha,cr,ce,mu_g,gamma]
        if not isinstance(parameters, list):
            raise TypeError('Input of myModel is of type list')

        if len(parameters) != 6:
            raise RuntimeError('Input list must be of length 6, containing [datsize, alpha, cr, ce, mu_g, gamma].')

        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector, name)

    def _check_input(self, input_values):
        # Check whether input has correct type or format
        if len(input_values) != 6:
            raise ValueError('Number of parameters of myModel must be 6.')

        # Check whether input is from correct domain (FIXME)
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]

        if not isinstance(datsize, np.int64):
            return False
        if datsize <= 0 or alpha < 0 or cr < 0 or ce < 0 or mu_g < 0 or gamma < 0:
            return False

        return True

    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output of the normal distribution is always a number.')
        return True

    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, k, rng=np.random.RandomState()):
        # Extract the input parameters
        datsize = input_values[0]
        alpha = input_values[1]
        cr = input_values[2]
        ce = input_values[3]
        mu_g = input_values[4]
        gamma = input_values[5]
        
        # Do the actual forward simulation
        vector_of_k_samples = self.mymodel_combined(datsize, alpha, cr, ce, mu_g, gamma, k) # For combined simulation
	# Format the output to obey API
        result = [np.array([x]) for x in vector_of_k_samples]
        return result

    

    def mymodel_combined(self, datsize, alpha, cr, ce, mu_g, gamma, k):

        # import warnings
        # warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
        
        weight = 1/np.sqrt(30)
        result = np.zeros((6,11))
        for ind_k in range(k):
            #Simulate at 3 Gy (clonogenic)
            fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,3.)
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            result[0,10] = np.log10(Ps[-1,0])
            
            # Simulate at 6 Gy (comet and clonogenic and comet)
            fun=lambda t, x: dP_dt(t,x,alpha,cr,ce,mu_g,gamma,6.)
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            Ps = odeint(fun, P0, ts, hmax = 1./3600 )
            result[1,10] = np.log10(Ps[-1,0])
            
            result[0,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes( 15) )
            result[1,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes( 30) )
            result[2,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes( 60) )
            result[3,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes(120) )
            result[4,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes(240) )
            result[5,0:10] = weight*get_normalized_row(Ps, get_index_from_minutes(360) )
                
            result = [result]
	
        return result


