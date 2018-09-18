# this code calculates MEX potential
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special 
import scipy.integrate as integrate
import scipy.interpolate as interpolate

# function to load halo path and return values in spherical coordinates
# with theta and phi having correct phase 

def LoadHalo(loadhalopath):

	mass,x,y,z,vx,vy,vz = np.loadtxt(loadhalopath,skiprows=1,unpack= True) 
	# need to calculate the radius each partcile is at 
	r = (x**2+y**2+z**2)**(0.5)  
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

	return(mass, r, theta_pos, phi_pos)

# r is radius of partciles [array] from LoadHalo, n is number of bins in r, logarithimcally
#max_r is the max radius you want the halo to subtend to over all time steps
def bin_r_log(r,n,max_r):
	n1 = n + 1
	# max radial bin                                   
	max_bin = np.log10(max_r)
	# make logarithimically spaced radial bins
	linear_space = np.linspace(0,max_bin,n1)
	log_bins = 10**(linear_space)
	bincenters = 0.5*(log_bins[1:]+log_bins[:-1])
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

	return(bincenters,r_trunc,log_bins)

# calculate the radial basis functions 
# listed in list as (0,0), (1,-1), (1,0),(1,1),(2,-2) .......
def basis_function(l_max, bincenters, log_bins, r_trunc, mass, theta_pos, phi_pos):
	# make list of l,m
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

	# compute radial bin widths, for calculation of basis functions
	deltar = []
	for i in range(len(log_bins)-1):
		width = (log_bins[i+1] - log_bins[i])
		deltar.append(width)
	# calculate r^2*deltar for basis function calculation
	r_sqaured_deltar = deltar*bincenters**2

	# calculate basis functions
	basis_functions = []
	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:

				if m == 0 :
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
						
					basis_functions.append(density)

				elif m < 0 :

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

					basis_functions.append(density)


				elif m > 0:

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

					basis_functions.append(density)


	return(basis_functions)

# calculate the potential now function 
# theta, phi and r_position is where you want the potential gven for 
# l_max must be the same as basis_function 
def MEX_potential_point(theta, phi , r_position, l_max, bincenters, basis_functions):
	G = 1

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

	# BT 2-122 1st ed. calculation

	phi_contribution = []

	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:


				if m < 0 :
				
				
					# function for integral1 
					func1 = basis_functions[count]*bincenters**(l+2)
					# function for integral2 
					func2 = basis_functions[count]/(bincenters**(l-1))


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

					# function for integral1 
					func1 = basis_functions[count]*bincenters**(l+2)
					# function for integral2 
					func2 = basis_functions[count]/(bincenters**(l-1))
					
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
					
					# function for integral1 
					func1 = basis_functions[count]*bincenters**(l+2)
					# function for integral2 
					func2 = basis_functions[count]/(bincenters**(l-1))

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

			
				count += 1

			else:

				zero_array = np.zeros(len(bincenters[1:]))
				phi_contribution.append(zero_array)

	# sum each contribtuion term by term to get total potential
	phi_contribution_sum = np.sum(phi_contribution, axis=0)

	## interpolate data to return potential at point wanted
	interpolate_potential = interpolate.interp1d(bincenters[1:],phi_contribution_sum, kind='cubic')

	return(interpolate_potential(r_position), phi_contribution_sum, phi_contribution)


def MEX_potential_plot(theta, phi, l_max, bincenters, basis_functions):
	G = 1

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

	# BT 2-122 1st ed. calculation

	phi_contribution = []

	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:


				if m < 0 :
				
				
					# function for integral1 
					func1 = basis_functions[count]*bincenters**(l+2)
					# function for integral2 
					func2 = basis_functions[count]/(bincenters**(l-1))


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

					# function for integral1 
					func1 = basis_functions[count]*bincenters**(l+2)
					# function for integral2 
					func2 = basis_functions[count]/(bincenters**(l-1))
					
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
					
					# function for integral1 
					func1 = basis_functions[count]*bincenters**(l+2)
					# function for integral2 
					func2 = basis_functions[count]/(bincenters**(l-1))

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

			
				count += 1

			else:

				zero_array = np.zeros(len(bincenters[1:]))
				phi_contribution.append(zero_array)

	# sum each contribtuion term by term to get total potential
	phi_contribution_sum = np.sum(phi_contribution, axis=0)

	# plot of potential from integral and from interpolation 																
	plt.plot(bincenters[1:],phi_contribution_sum)								
	plt.xlabel(r"$r$ [kpc]")											
	plt.ylabel(r'$\Phi$')	
	plt.title(r"MEX Potential, $l=%s$, $\theta = %s$, $\phi = %s$"  %(l_max,theta,phi))
	#plt.savefig(savefilepath)
	plt.show()
	plt.clf()




				

























