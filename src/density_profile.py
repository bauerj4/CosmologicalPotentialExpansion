# G=1 

# this code will calculate the density profile of galaxies 
# will also calctlate analytic potential and density 


# importing useful packages 
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scip
plt.rcParams['agg.path.chunksize'] = 100000000
import scipy.special as special 


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

# number of bins for radial binning
n = 28
#number of bins used to find fit paramters for analytic profile
k = 28

# file paths for saving density profile plots
#chnage file paths as required  

#log(rho*r^2) vs log(r) format
savefilepath1 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s_format1.eps" %(galaxy,galaxy)
#log(rho) vs log(r) format
savefilepath2 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s_format2.eps" %(galaxy,galaxy)
# rho vs r format
savefilepath3 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s_format3.eps" %(galaxy,galaxy)
# best fit of desity profile overlayed log(rho) vs log(r) format ignoring points
savefilepath4 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s_densityprofilefitt_pointsecluded.eps" %(galaxy,galaxy)
# best fit of desity profile overlayed log(rho) vs log(r) format without ignoring points
savefilepath5 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s_densityprofilefitt_nopointsexcluded.eps" %(galaxy,galaxy)
# integral of potential spherical case
savefilepath6 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/density_profile_%s_potential.eps" %(galaxy,galaxy)

# MEX potentil (Binney and Tremaine integral (2-122), 1st Ed.)

# format for sph_harm ... sph_harm(m,l,theta,phi)
# theta runs from 0 to 2pi
# phi runs from 0 to pi 

G = 1
l_max = 0

# position on galaxy (in radians)
theta = 0
phi = 0



########################################################################################################################################################################

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

#compute the density in each radial shell 
density = counts/volumes 
#density is in units of #particles/kpc^3 

# plot the density profiles 

# log(rho*r^2) vs log(r) 
plt.plot(bincenters,density*bincenters**2,"o")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$\log(r)$")
plt.ylabel("$log(r^2$"r'$\rho$'")")
plt.title("Density profile of %s" %(galaxy))
plt.savefig(savefilepath1)
#plt.show()
plt.clf()

# log(rho) vs log(r)
plt.plot(bincenters,density,"o")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("$\log(r)$")
plt.ylabel("$log($"r'$\rho$'")")
plt.title("Density profile of %s" %(galaxy))
plt.savefig(savefilepath2)
#plt.show()
plt.clf()

# rho vs r format
plt.plot(bincenters,density,"o")
plt.xlabel("$r$")
plt.ylabel(r'$\rho$')
plt.title("Density profile of %s" %(galaxy))
plt.savefig(savefilepath3)
#plt.show()
plt.clf()

# fitting density profile to find rho_0 and R_s function

def func(x,rho_0,R_s):
	Rho = (rho_0/((x/R_s)*(1+(x/R_s))**2))
	return Rho 


# fit ignoring k points at end 
popt1, pcov1 = scip.curve_fit(func, bincenters[0:k], density[0:k])



# produce analaytic density profile overlayed binned one

# parameters from NFW profile pulled out of optimization
rho_01 = popt1[0]
R_s1 = popt1[1]

# arrange r in ascending order 
r_sort = sorted(r)
Rho1 = (rho_01/((r_sort/R_s1)*(1+(r_sort/R_s1))**2))

# plot of log(r) vs log(rho) profile with fit 
plt.plot(r_sort,Rho1,label=r'NFW profile, $R_s$= %5.5f,'r'$\rho$''$_0$= %5.5f' %(R_s1,rho_01) )
plt.plot(bincenters[0:k],density[0:k],"o",label='Binned Data')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$\log(r)$")
plt.ylabel(r"$\log($"r'$\rho$'")")
plt.legend()
plt.title("Density profile %s with fit, points ignored" %(galaxy))
plt.savefig(savefilepath4)
#plt.show()
plt.clf()

# fit without ignoring points 

popt2, pcov2 = scip.curve_fit(func, bincenters, density)

rho_02 = popt2[0]
R_s2 = popt2[1]

Rho2 = (rho_02/((r_sort/R_s2)*(1+(r_sort/R_s2))**2))

# plot of log(r) vs log(rho) profile with fit 
plt.plot(r_sort,Rho2,label=r'NFW profile, $R_s$= %5.5f,'r'$\rho$''$_0$= %5.5f' %(R_s2,rho_02) )
plt.plot(bincenters,density,"o",label='Binned Data')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$\log(r)$")
plt.ylabel(r"$\log($"r'$\rho$'")")
plt.legend()
plt.title("Density profile %s with fit, no points ignored" %(galaxy))
plt.savefig(savefilepath5)
#plt.show()
plt.clf()

#### MEX potential (specific case l=0) ####

# function for integral1 
func1 = density*bincenters**2
# function for integral2 
func2 = density*bincenters

#numerically integrating each term using trapezoid rule

# not needed trapcenters *********************
#need to find center of trapezoids 
trapcenters = 0.5*(bincenters[1:]+bincenters[:-1])

# Term 1 of potenetil at each bincenter
x = 0 
int1vals = []
int1vals_running_sum = []
while x != len(bincenters)-1 : 
 	value1 = ((func1[x+1]+func1[x])/2)*(bincenters[x+1]-bincenters[x])
 	int1vals.append(value1)
 	sum_int1vals = sum(int1vals)
 	int1vals_running_sum.append(sum_int1vals)
 	x += 1 


# Term 2 of potential at each bincenter
x = 0 
int2vals = []
int2vals_running_sum = []
while x != len(bincenters)-1 : 
	value2 = ((func1[x+1]+func1[x])/2)*(bincenters[x+1]-bincenters[x])
	int2vals.append(value1)
	sum_int2vals = sum(int2vals)
	int2vals_running_sum.append(sum_int2vals)
	x += 1 

# now we compute the value of the integrals combining term 1 and 2 
x = 0 
combinedvals = []
while x != len(bincenters)-1 : 
	combval = int1vals_running_sum[x]/bincenters[x+1] + int2vals_running_sum[-x-1]
	combinedvals.append(combval)
	x +=1

combinedvals = np.array(combinedvals)
Phi_r = -4*np.pi*G*(combinedvals)*np.real(special.sph_harm(0,0,theta,phi))



#calculate analytic potential now using NFW profile fit parameters #
                                                                    #
def potential(x,rho_0,R_s):                                        #
	G=1																#
	phi_analytic = -(4*np.pi*G*rho_0*R_s**3*np.log(1+x/R_s))/x		#		 
	return phi_analytic												#
																	#
phi_analytic = potential(r,rho_02,R_s2)							#
																	#

# plot of potential from integral and from fit (not saved)																	#
plt.plot(r,phi_analytic,"o",label='from fit')
plt.plot(bincenters[1:],Phi_r,label='from integral')										#
plt.xlabel(r"$r$ [kpc]")											#
plt.ylabel(r'$\Phi$')	
plt.title("Analytic potential and from integral, $l=0$")
plt.legend(fontsize = '8')
plt.show()
plt.clf()
											

#plot of potential from integral
plt.plot(bincenters[1:],Phi_r)
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r'$\Phi$')
plt.title("from integral, $l=0$ case")
plt.savefig(savefilepath6)
#plt.show()
plt.clf()

# general case (l not neccicarily zero)


































# old routines not currently needed 


#####################################################################
# calculate analytic potential now using NFW profile fit parameters #
                                                                    #
#def potential(x,rho_0,R_s):                                        #
#	G=1																#
#	phi_analytic = -(4*np.pi*G*rho_0*R_s**3*np.log(1+x/R_s))		#		 
#	return phi_analytic												#
																	#
#phi_analytic = potential(r,rho_01,R_s1)							#
																	#
																	#
#plt.plot(r,phi_analytic,"o")										#
#plt.xlabel(r"$r$ [kpc]")											#
#plt.ylabel(r'$\Phi$')												#
#plt.title("from fit of density paramaters")						#
#plt.show()															#
#plt.clf()															#
#####################################################################


###################################
# Check of bin centers and counts #
# not included in code normally   #
#plt.plot(bincenters,counts)      #
#plt.hist(r,bins=log_bins)        #
#plt.show()                       #
#plt.clf()                        #
###################################



