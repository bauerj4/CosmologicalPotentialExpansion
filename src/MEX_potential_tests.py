# G=1 

# this code will calculate MEX potential 

# importing useful packages 
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scip
plt.rcParams['agg.path.chunksize'] = 100000000
import scipy.special as special 
import scipy.integrate as integrate

# controls the parameters of the density profile edited for specific halos
# change as needed

# importing halo initial conditions 
loadhalopath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/halos/model_A_spherical_halo_BW2018b.ascii"
m,x,y,z,vx,vy,vz = np.loadtxt(loadhalopath,skiprows=1,unpack= True) 
# m is mass in units of 2.302e9 solar masses
# x,y,z is position in kpc 
# vx,vy,vz is velocity in units of 100 km/s


# name of galaxy simulation for title of plots / savefilepath
galaxy = "model_A_spherical_halo_BW2018b"

# number of bins for radial binning ( # of bins is n-1)
n = 31


# file paths for saving density profile fit and parmaters / potnential
#chnage file paths as required  


# MEX potentil (Binney and Tremaine integral (2-122), 1st Ed.)

# format for sph_harm ... sph_harm(m,l,theta,phi)
# theta runs from 0 to 2pi
# phi runs from 0 to pi 

G = 1
l_max = 2

# position on galaxy (in radians)
theta = 0
phi = 0


########################################################################################################################################################################

##### need to radially bin the the particles logarithimically #####

# need to calculate the radius each partcile is at 
r = (x**2+y**2+z**2)**(0.5)                        
# r is the radius of the particles in units of kpc 
                                                   
# max radial bin                                   
max_bin = np.log10(max(r))                               

# make logarithimically spaced bins for histograms 
linear_space = np.linspace(0,max_bin,n)
log_bins = 10**(linear_space)
bincenters = 0.5*(log_bins[1:]+log_bins[:-1])

# bin the data
counts, bins, patches = plt.hist(r,bins=log_bins)
plt.clf()
# counts is the counts in each bin


##############################################################
# compute volume conained in each radial bin 
volumes = []
#volumes is volume contained in each radial bin [kpc^3]
for i in range(len(log_bins)-1):
	V = (4/3)*np.pi*(log_bins[i+1]**3 - log_bins[i]**3)
	volumes.append(V)

#compute the density in each radial shell 
density = counts/volumes 
#density is in units of #particles/kpc^3 
#############################################################

















