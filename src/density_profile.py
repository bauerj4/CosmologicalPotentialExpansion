# G=1 

# this code will calculate the density profile of galaxies 

# importing useful packages 
import numpy as np
import matplotlib.pyplot as plt



# controls the parameters of the density profile edited for specific halos
# change as needed

# importing halo initial conditions 
loadhalopath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/halos/model_A_spherical_halo_BW2018b.ascii"
m,x,y,z,vx,vy,vz = np.loadtxt(loadhalopath,skiprows=1,unpack= True) 
# m is mass in units of 2.302e9 solar masses
# x,y,z is position in kpc 
# vx,vy,vz is velocity in units of 100 km/s


# name of galaxy simulation for title of density profile / savefilepath
galaxy = "model_A_spherical_halo_BW2018b"

# number of bins for radial binning
n = 30 

# file path for saving density profile 
savefilepath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s.eps" %(galaxy,galaxy)





#################### NO USER EDITING REQUIRED PAST THIS POINT ##############

##### need to radially bin the the particles logarithimically #####

# need to calculate the radius each partcile is at 
r = (x**2+y**2+z**2)**(0.5)                        
# r is the radius of the particles in units of kpc 
                                                   
# max radial bin                                   
m = np.log10(max(r))                               

# make logarithimically spaced bins for histograms 
linear_space = np.linspace(0,m,n)
log_bins = 10**(linear_space)
bincenters = 0.5*(log_bins[1:]+log_bins[:-1])

# bin the data 
counts, bins, patches = plt.hist(r,bins=log_bins)
plt.clf()
# counts is the counts in each bin

# compute volume conained in each radial bin 
volumes = []
#volumes is volume contained in each radial bin [kpc^3]
for i in range(len(log_bins)-1):
	V = (4/3)*np.pi*(log_bins[i+1]**3 - log_bins[i]**3)
	volumes.append(V)

#compute the density in each radia shell 
density = counts/volumes 
#density is in units of #particles/kpc^3 

# plot the density profile 
plt.plot(bincenters,density*bincenters**2,"o")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$\log(r)$")
plt.ylabel("$log(r^2$"r'$\rho$'")")
plt.title("Density profile of %s" %(galaxy))
plt.savefig(savefilepath)
#plt.show()
plt.clf()



# Check of bin centers and counts #
# not included in code normally   #
#plt.plot(bincenters,counts)      #
#plt.hist(r,bins=log_bins)        #
#plt.show()                       #
#plt.clf()                        #
###################################



