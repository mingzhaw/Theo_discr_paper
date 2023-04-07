import numpy as np

############################################
#                                          #
#   Clusters of hits with n_target = 200   #
#                                          #
############################################


n_target = int(200)
q = 0.9

n_steps = 12 #Steps of 0.5 Gy till 6 Gy

N_cells = int(1e6) #Number of cells simulated
population = np.zeros(N_cells) #Vector to save the final cell population

#Definition of vectors to study the distribution of cluster with diverse hits and the obtained synthetic comet
n_hits = [1, 2, 3, 4, 5, 6, 7, 8, 9]
cluster_distribution = np.zeros(len(n_hits))
population_freq = np.zeros(len(n_hits))


#Loop to simulate the generation of clusters
for i in range(0,N_cells):
    #Vector to save the accumulation of hits in the different target sites
    hits = np.zeros(n_target)
    
    #Loop to generate damage within the cell
    for j in range(0,n_steps):
        #Random number generation to decide if a hit is generaed in each target site
        r = np.random.uniform(0,1,n_target)
        hits[r>=q] += 1
    
    #Save the cluster with the highest number of hits (i.e., the population to which the cell belongs) and the contribution to the final synthetic comet                        
    population[i] = np.max(hits)
    population_freq[int(np.max(hits)-1)] += 1/N_cells
    
    #Calculate the hit distributions
    for j in range(0,len(cluster_distribution)):
        cluster_distribution[j] += len(hits[hits == n_hits[j]])/N_cells
        

#Save the results
np.savetxt('clusters_200.dat',np.transpose([n_hits, cluster_distribution, population_freq]), header="cluster amount population", comments='')




###############################################
#                                             #
#   Clusters of hits with n_target = 10.000   #
#                                             #
###############################################

n_target = int(1e4)
q = 0.998

n_steps = 12 #Steps of 0.5 Gy till 6 Gy

N_cells = int(1e6) #Number of cells simulated
population = np.zeros(N_cells) #Vector to save the final cell population

#Definition of vectors to study the distribution of cluster with diverse hits and the obtained synthetic comet
n_hits = [1, 2, 3, 4, 5, 6, 7, 8, 9]
cluster_distribution = np.zeros(len(n_hits))
population_freq = np.zeros(len(n_hits))


#Loop to simulate the generation of clusters
for i in range(0,N_cells):
    #Vector to save the accumulation of hits in the different target sites
    hits = np.zeros(n_target)
    
    #Loop to generate damage within the cell
    for j in range(0,n_steps):
        #Random number generation to decide if a hit is generaed in each target site
        r = np.random.uniform(0,1,n_target)
        hits[r>=q] += 1
    
    #Save the cluster with the highest number of hits (i.e., the population to which the cell belongs) and the contribution to the final synthetic comet                        
    population[i] = np.max(hits)
    population_freq[int(np.max(hits)-1)] += 1/N_cells
    
    #Calculate the hit distributions
    for j in range(0,len(cluster_distribution)):
        cluster_distribution[j] += len(hits[hits == n_hits[j]])/N_cells



#Save the results
np.savetxt('clusters_10_4.dat',np.transpose([n_hits, cluster_distribution, population_freq]), header="cluster amount population", comments='')
