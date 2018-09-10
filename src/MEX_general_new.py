# once written into function call
# user inputs # of bins, (radius, theta, phi) the position they want the potential 
# to be retunred in spherical coordiinates and force to be retunred in cartesian
# user also inputs the load halo path
# user also tells lmax

# G=1 

# this code will calculate MEX potential 

# importing useful packages 
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scip
plt.rcParams['agg.path.chunksize'] = 100000000
import scipy.special as special 
import scipy.integrate as integrate
import scipy.interpolate as interpolate
# controls the parameters of the density profile edited for specific halos
# change as needed

# importing halo initial conditions 
loadhalopath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/halos/model_A_spherical_halo_BW2018b.ascii"
mass,x,y,z,vx,vy,vz = np.loadtxt(loadhalopath,skiprows=1,unpack= True) 
# m is mass in units of 2.302e9 solar masses
# x,y,z is position in kpc 
# vx,vy,vz is velocity in units of 100 km/s

# name of galaxy simulation for title of plots / savefilepath
galaxy = "model_A_spherical_halo_BW2018b"

# number of bins for radial binning ( # of bins is n-1) (logarithmic)
n1 = 33

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
l_max = 5
theta = 0
phi = 0

# position on galaxy for potential is in radians
# the potential will be dined along the directions of the theta and phi bin centers
#theta and phi bins must have same dimension, as well as same dimesnion of radial bins


########################################################################################################################################################################

##### need to radially bin the the particles logarithimically #####

# ordering of points is preserved becasue the ordering is not touched

# need to calculate the radius each partcile is at 
r = (x**2+y**2+z**2)**(0.5)                        
# r is the radius of the particles in units of kpc 
                                                   
# max radial bin                                   
max_bin = np.log10(500)  #np.log10(max(r))                               

# make logarithimically spaced radial bins
linear_space = np.linspace(0,max_bin,n1)
log_bins = 10**(linear_space)
bincenters = 0.5*(log_bins[1:]+log_bins[:-1])

# compute radial bin widths, for calculation of basis functions
deltar = []
for i in range(len(log_bins)-1):
	width = (log_bins[i+1] - log_bins[i])
	deltar.append(width)
# calculate r^2*deltar for basis function calculation
r_sqaured_deltar = deltar*bincenters**2

# calculate positions of particle in spherical coordinates (theta_pos and phi_pos)
theta_pos = np.arctan2(y,x)

phi_pos = np.arccos(z/r)

# covert theta positions from -pi to pi --> 0 to 2pi 

theta_pos_correct_phase = []

for x in theta_pos:
	if x < 0:
		X = x + 2*np.pi 
		theta_pos_correct_phase.append(X)
	elif x > 0:
		theta_pos_correct_phase.append(x)
	elif x == 0:
		theta_pos_correct_phase.append(x)

theta_pos = theta_pos_correct_phase

# bin the data logarithically in r 
# coordinates for binning go as (r,phi,theta) (math convention not physics, due to scipy definitions)
# radial bins made above 
# convert the r data to the corresponding bin center the ppoint falls in, only for r

r_trunc = r
i = 0
while i < len(bincenters):

	count = 0
	for x in r: 

		if log_bins[i] <= x < log_bins[i+1]:

			r_trunc[count] = bincenters[i]


		else:
			pass 
		
		count += 1 		

	i += 1



#calculate generally the potential

# make list of l
a = 0 
L = []
while a <= l_max:
	L.append(a)
	a += 1
L = np.array(L)

# make list of m 
b = -l_max 
M = []
while b <= l_max:
	M.append(b)
	b += 1
M = np.array(M)

# compute equation (Binney and Tremaine integral (2-122), 1st Ed.), as well radial basis functions
# defined in latex document 

#phi_contribution is the contribution to potnetial from each l,m possibility
# condition on m determines if statements due to spherical harmonics 




# for testing
#L = [1]
#M = [0]
##############



phi_contribution = []

for l in L:
	for m in M: 
		

		if np.absolute(m) <= l:

			if m < 0 :
				
				# calculate radial basis function for l,m

				density = []
				i = 0 
				while i < len(bincenters):

					density_cont = []
	
					count = 0
					for x in r_trunc:

						if x == bincenters[i]:

							term = (1/r_sqaured_deltar[i])*mass[count]* np.real(((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta_pos[count],phi_pos[count])-(-1)**m*special.sph_harm(m,l,theta_pos[count],phi_pos[count]))).conjugate())
							density_cont.append(term)



						else:
							pass

						count += 1 

					density.append(sum(density_cont))

					i += 1 
				
				
				# function for integral1 
				func1 = density*bincenters**(l+2)
				# function for integral2 
				func2 = density/(bincenters**(l-1))


				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 needs to be calculated in following loop
				

				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals = []
				while x != len(bincenters)-1 : 
					int2valrunningsum = integrate.cumtrapz(func2[x:],bincenters[x:])
					combval = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[len(int2valrunningsum)-1]*bincenters[x+1]**(l)
					combinedvals.append(combval)
					x +=1 
				combinedvals = np.array(combinedvals)

				#compute the potential term for this l and m comb
				phi_comb = (-4*np.pi*G*combinedvals*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)
				# append all the contributions of phi into another list of the phi
				phi_contribution.append(phi_comb)


			

			elif m == 0 :
				
				# calculate radial basis function for l,m
				
				density = []
				i = 0 
				while i < len(bincenters):

					density_cont = []
	
					count = 0
					for x in r_trunc:

						if x == bincenters[i]:

							term = (1/r_sqaured_deltar[i])*mass[count]*np.real((special.sph_harm(m,l,theta_pos[count],phi_pos[count])).conjugate())

							density_cont.append(term)



						else:
							pass

						count += 1 

					density.append(sum(density_cont))

					i += 1 
				

				# function for integral1 
				func1 = density*bincenters**(l+2)
				# function for integral2 
				func2 = density/(bincenters**(l-1))
				
				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 needs to be calculated in following loop
		

				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals = []

				while x != len(bincenters)-1 : 
					int2valrunningsum = integrate.cumtrapz(func2[x:],bincenters[x:])
					combval = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[len(int2valrunningsum)-1]*bincenters[x+1]**(l)
					combinedvals.append(combval)
					x +=1 
				combinedvals = np.array(combinedvals)
				#compute the potential term for this l and m comb
				phi_comb =(-4*np.pi*G*combinedvals*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)
				# append all the contributions of phi into another list of the phi
				phi_contribution.append(phi_comb)


			elif m > 0:
				
				# calculate radial basis function for l,m

				density = []
				i = 0 
				while i < len(bincenters):

					density_cont = []
	
					count = 0
					for x in r_trunc:

						if x == bincenters[i]:

							term = (1/r_sqaured_deltar[i])*mass[count]* np.real(((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta_pos[count],phi_pos[count])+(-1)**m*special.sph_harm(-m,l,theta_pos[count],phi_pos[count]))).conjugate())

							density_cont.append(term)



						else:
							pass

						count += 1 

					density.append(sum(density_cont))

					i += 1 
				


				# function for integral1 
				func1 = density*bincenters**(l+2)
				# function for integral2 
				func2 = density/(bincenters**(l-1))

				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 needs to be calculated in following loop
				
				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals = []
				while x != len(bincenters)-1 : 
					int2valrunningsum = integrate.cumtrapz(func2[x:],bincenters[x:])
					combval = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[len(int2valrunningsum)-1]*bincenters[x+1]**(l)
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

print(int2valrunningsum)

# sum each contribtuion term by term to get total potential
phi_contribution_sum = np.sum(phi_contribution, axis=0)

# interpolate data to return potential at point wanted
interpolate_potential = interpolate.interp1d(bincenters[1:],phi_contribution_sum, kind='cubic')

#potential_returned = interpolate_potential(radius)


xnew = np.linspace(bincenters[1],max(r_trunc),100)

#print(phi_contribution_sum)



# plot of potential from integral and from interpolation 																
plt.plot(bincenters[1:],phi_contribution_sum)
plt.plot(xnew,interpolate_potential(xnew))										
plt.xlabel(r"$r$ [kpc]")											
plt.ylabel(r'$\Phi$')	
plt.title("MEX, $l=0$")
plt.xlim(0,max(r_trunc))
plt.show()
plt.clf()




























# loop to calculate the density function  (will need to be inserted into other loop(s))

#density = []
#i = 0 
#while i < len(bincenters):
#
#	density_cont = []
#	
#	count = 0
#	for x in r_trunc:
#
#		if x == bincenters[i]:
#
#			term = (1/r_sqaured_deltar[i])*m[count]* (spherical harmonics depending on condition)
#
#			density_cont.append(term)
#
#
#		else:
#			pass
#
#		count += 1 
#	density.append(sum(density_cont))
#
#	i += 1 


