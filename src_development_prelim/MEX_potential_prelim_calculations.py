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

print(max(m))
print(min(m))
print(m)
# name of galaxy simulation for title of plots / savefilepath
galaxy = "model_A_spherical_halo_BW2018b"

# number of bins for radial binning ( # of bins is n-1)
n = 31


# file paths for saving density profile fit and parmaters / potnential
#chnage file paths as required  

# best fit of desity profile overlayed log(rho) vs log(r) format without ignoring points
savefilepath1 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/MEX_potential_code_plots/densityprofilefit.eps" %(galaxy)
# integral of potential spherical case
savefilepath2 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/MEX_potential_code_plots/potential_int.eps" %(galaxy)
# potneital from integral anf from denisty profile fit
savefilepath3 = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/%s_plots/MEX_potential_code_plots/potential_int_analytic_overlayed.eps" %(galaxy)


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

# compute volume conained in each radial bin 
volumes = []
#volumes is volume contained in each radial bin [kpc^3]
for i in range(len(log_bins)-1):
	V = (4/3)*np.pi*(log_bins[i+1]**3 - log_bins[i]**3)
	volumes.append(V)

#compute the density in each radial shell 
density = counts/volumes 
#density is in units of #particles/kpc^3 


# fitting density profile to find rho_0 and R_s function

def func(x,rho_0,R_s):
	Rho = (rho_0/((x/R_s)*(1+(x/R_s))**2))
	return Rho 


popt2, pcov2 = scip.curve_fit(func, bincenters, density)

rho_02 = popt2[0]
R_s2 = popt2[1]

# arrange r in ascending order 
r_sort = sorted(r)

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
#plt.savefig(savefilepath1)
#plt.show()
plt.clf()



#### MEX potential (specific case l=0)

# function for integral1 
func1 = density*bincenters**2
# function for integral2 
func2 = density*bincenters

#numerically integrating each term using trapezoid rule

#term 1

int1valrunningsum = integrate.cumtrapz(func1,bincenters)

# term 2 

int2valrunningsum = integrate.cumtrapz(func2,bincenters)

# now we compute the value of the integrals combining term 1 and 2 
x = 0 
combinedvals = []
while x != len(bincenters)-1 : 
	combval = int1valrunningsum[x]/bincenters[x+1] + int2valrunningsum[-x-1]
	combinedvals.append(combval)
	x +=1

combinedvals = np.array(combinedvals)
Phi_r = -4*np.pi*G*(combinedvals)*np.real(special.sph_harm(0,0,theta,phi))

						

#calculate analytic potential now using NFW profile fit parameters  
                                                                    
def potential(x,rho_0,R_s):                                         
	G=1																
	phi_analytic = -(4*np.pi*G*np.real(special.sph_harm(0,0,theta,phi))*rho_0*R_s**3*np.log(1+x/R_s))/x				 
	return phi_analytic											
																	
phi_analytic = potential(r,rho_02,R_s2)							    
																	

# plot of potential from integral and from fit																
plt.plot(r,phi_analytic,'o',label='from fit')
plt.plot(bincenters[1:],Phi_r,label='from integral')										
plt.xlabel(r"$r$ [kpc]")											
plt.ylabel(r'$\Phi$')	
plt.title("Analytic potential and from integral, $l=0$")
plt.legend(fontsize = '8')
#plt.savefig(savefilepath3)
plt.show()
plt.clf()
											

#plot of potential from integral
plt.plot(bincenters[1:],Phi_r)
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r'$\Phi$')
plt.title("from integral, $l=0$ case")
#plt.savefig(savefilepath2)
plt.show()
plt.clf()


#calculate generally the potential

# make list of l
a = 0 
L = []
while a <= l_max:
	L.append(a)
	a += 1
L = np.array(L)

# need to check limit of m in real sperical harmonic case

# make list of m 
b = -l_max 
M = []
while b <= l_max:
	M.append(b)
	b += 1
M = np.array(M)

######### just for testing
L= [1]
M= [0]

##########

# compute equation (Binney and Tremaine integral (2-122), 1st Ed.)

#phi_contribution is the contribution to potnetial from each l,m possibility
phi_contribution = []

for l in L:
	for m in M: 
		
		# function for integral1 
		func1 = density*bincenters**(l+2)
		# function for integral2 
		func2 = density/bincenters**(l-1)

		if np.absolute(m) <= l:

			if m < 0 :
		
				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 
				int2valrunningsum = integrate.cumtrapz(func2,bincenters)

				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals = []
				while x != len(bincenters)-1 : 
					combval = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[-x-1]*bincenters[x+1]**(l)
					combinedvals.append(combval)
					x +=1 
				combinedvals = np.array(combinedvals)

				#compute the potential term for this l and m comb
				phi_comb = (-4*np.pi*G*combinedvals*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)
				# append all the contributions of phi into another list of the phi
				phi_contribution.append(phi_comb)


			

			elif m == 0 :

				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 
				int2valrunningsum = integrate.cumtrapz(func2,bincenters)

				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals = []

				while x != len(bincenters)-1 : 
					combval = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[-x-1]#*bincenters[x+1]**(l)
					combinedvals.append(combval)
					x +=1 
				combinedvals = np.array(combinedvals)
				#compute the potential term for this l and m comb
				phi_comb =(-4*np.pi*G*combinedvals*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)
				# append all the contributions of phi into another list of the phi
				phi_contribution.append(phi_comb)


			elif m > 0:

				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 
				int2valrunningsum = integrate.cumtrapz(func2,bincenters)

				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals = []
				while x != len(bincenters)-1 : 
					combval = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[-x-1]#*bincenters[x+1]**(l)
					combinedvals.append(combval)
					x +=1 
				combinedvals = np.array(combinedvals)
				#compute the potential term for this l and m comb
				phi_comb = (-4*np.pi*G*combinedvals*np.real((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)
				# append all the contributions of phi into another list of the phi
				phi_contribution.append(phi_comb)

		else:

			zero_array = np.zeros(len(bincenters[1:]))
			phi_contribution.append(zero_array)


phi_contribution_sum = np.sum(phi_contribution, axis=0)






plt.plot(bincenters[1:],phi_contribution_sum)
plt.plot(bincenters[1:],Phi_r,'o')
plt.xlabel(r"$r$ [kpc]")
plt.ylabel(r'$\Phi$')
plt.title("from integral, $l=1$ case")
#plt.savefig(savefilepath2)
#plt.show()
plt.clf()


#print(phi_contribution_sum)
#print("")
#print(Phi_r)
#print(L)
#print(M)









#print(phi_contribution)	
#print(len(phi_contribution))
#print(Phi_r)





















































#check of potential integral, laplacian to obtain back potential 
# from origianl PDE of MEX expansion
#a = np.diff(Phi_r,n=1)
#b = a*bincenters[2:]**2
#c = np.diff(b,n=1)
#d = c/(bincenters[3:]**(2))
#rho_sanitycheck = d/(4*np.pi*G)

#print(len(bincenters[1:]))
#print(len(c))
#print(a)

#plt.plot(bincenters[3:],rho_sanitycheck)
#plt.plot(bincenters,density)
#plt.xlim(60,140)
#plt.show()



















