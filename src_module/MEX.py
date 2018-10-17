# this code calculates MEX potential
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special 
import scipy.integrate as integrate
import scipy.interpolate as interpolate




# function to load halo path and return values in spherical coordinates
# with theta and phi having correct phase 
 
# need to put in fudge factors on theta and phi positions for now *************

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


# as each position has the same int_lm, we will calculate int_lm and then sum them to calculate the potential
# this int is for potential and theta and phi force directions 
def MEX_integral_lm_potential(l_max,basis_functions,bincenters):

	
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

	# calculate each integral contribution for lm, to be summed in potential function 
	int_lm = []

	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:
				
				
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

				# append combinevals into int_lm

				int_lm.append(combinedvals)


				count += 1

			else:

				pass

	return(int_lm)


# this int_lm only works for the force in the radial direction only
def MEX_integral_lm_force_r(l_max,basis_functions,bincenters):

	
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

	# calculate each integral contribution for lm, to be summed in potential function 
	int_lm = []

	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:
				
				
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
					combval = -(l+1)*(int1valrunningsum[x]/bincenters[x+1]**(l+2)) + l*int2valrunningsum[len(int2valrunningsum)-1]*bincenters[x+1]**(l-1)
					combinedvals.append(combval)
					x +=1 
				combinedvals = np.array(combinedvals)

				# append combinevals into int_lm

				int_lm.append(combinedvals)


				count += 1

			else:

				pass

	return(int_lm)


# returns inetgarls for both force and potentials , int_lm1 is for potential and theta/phi forces ... int_lm2 is for radial force only
def MEX_integral_lm_force_potential(l_max,basis_functions,bincenters):

	
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

	# calculate each integral contribution for lm, to be summed in potential function 
	int_lm1 = []
	int_lm2 = []

	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:
				
				
				# function for integral1 
				func1 = basis_functions[count]*bincenters**(l+2)
				# function for integral2 
				func2 = basis_functions[count]/(bincenters**(l-1))

				#term 1
				int1valrunningsum = integrate.cumtrapz(func1,bincenters)

				# term 2 needs to be calculated in following loop
					

				# now we compute the value of the integrals combining term 1 and 2
				x = 0 
				combinedvals1 = []
				combinedvals2 = []
				while x != len(bincenters)-1 : 
					int2valrunningsum = integrate.cumtrapz(func2[x:],bincenters[x:])
					combval1 = int1valrunningsum[x]/bincenters[x+1]**(l+1) + int2valrunningsum[len(int2valrunningsum)-1]*bincenters[x+1]**(l)
					combval2 = -(l+1)*(int1valrunningsum[x]/bincenters[x+1]**(l+2)) + l*int2valrunningsum[len(int2valrunningsum)-1]*bincenters[x+1]**(l-1)
					combinedvals1.append(combval1)
					combinedvals2.append(combval2)
					x +=1 
				

				combinedvals1 = np.array(combinedvals1)
				combinedvals2 = np.array(combinedvals2)


				# append combinevals into int_lm

				int_lm1.append(combinedvals1)
				int_lm2.append(combinedvals2)


				count += 1

			else:

				pass

	return(int_lm1,int_lm2)





# calculate the potential now function 
# theta, phi and r_position is where you want the potential gven for 
# l_max must be the same as basis_function and MEX_integral_lm 
# (r_position, theta, phi) correspond to particle of interest 

def MEX_potential(l_max, r_position, theta, phi, bincenters, int_lm):
	
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
				

					#compute the potential term for this l and m comb
					phi_comb = (-4*np.pi*G*int_lm[count]*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)
					# append all the contributions of phi into another list of the phi
					phi_contribution.append(phi_comb)



				elif m == 0 :

					
					#compute the potential term for this l and m comb
					phi_comb =(-4*np.pi*G*int_lm[count]*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)
					# append all the contributions of phi into another list of the phi
					phi_contribution.append(phi_comb)


				elif m > 0:
					
					
					#compute the potential term for this l and m comb
					phi_comb = (-4*np.pi*G*int_lm[count]*np.real((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)
					# append all the contributions of phi into another list of the phi
					phi_contribution.append(phi_comb)

			
				count += 1

			else:

				pass

	# sum each contribtuion term by term to get total potential
	phi_contribution_sum = np.sum(phi_contribution, axis=0)

	## interpolate data to return potential at point wanted
	interpolate_potential = interpolate.interp1d(bincenters[1:],phi_contribution_sum, kind='cubic')
	# potential at particles position 
	PHI = interpolate_potential(r_position)

	return(PHI, phi_contribution_sum, phi_contribution)


# NEED TO FINISH ( correct factors on forces) ******


# calculate the forces on the particles 
# int order is same as the function that calculates it 
# (r_position, theta, phi, mass) correspond to the particle of interest
# int_lm1 is for potential and theta/phi forces ... int_lm2 is for radial force only
def MEX_force(l_max, r_position, theta, phi, mass, bincenters , int_lm1, int_lm2):

	# make function spherical harmonic derivative with respect to phi  (theta physics convention)
	def sph_harm_deriv(m,l,theta,phi):
		
		if l = m : 
			
			d_phi = m*(1/np.tan(phi))*special.sph_harm(m,l,theta,phi) 

		else:

			d_phi = m*(1/np.tan(phi))*special.sph_harm(m,l,theta,phi) + np.sqrt((l-m)*(l+m+1))*np.e**(-1j*theta)*special.sph_harm(m+1,l,theta,phi)

		return(d_phi)
	

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

	F_r = []
	F_theta = []
	F_phi = []


	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:


				if m < 0 :
				

					#compute the radial force
					f_r = (mass*4*np.pi*G*int_lm2[count]*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm1[count]*np.real(1j*m*(1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)

					#compute the phi force (theta in physics convention)
					f_phi = (4*mass*np.pi*G*int_lm1[count]*np.real((1/(1j*np.sqrt(2)))*(sph_harm_deriv(-m,l,theta,phi)-(-1)**m*sph_harm_deriv(m,l,theta,phi))))/(2*l+1)

					

					# append all the contributions of forces into another list of the forces
					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)



				elif m == 0 :

					
					#compute the radial force
					f_r =(4*mass*np.pi*G*int_lm2[count]*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm1[count]*np.real(1j*m*special.sph_harm(m,l,theta,phi)))/(2*l+1)

					#compute the phi force (theta in physics convention)
					f_phi = (4*mass*np.pi*G*int_lm1[count]*np.real(sph_harm_deriv(m,l,theta,phi)))/(2*l+1)


					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)



				elif m > 0:
					
					
					#compute the radial force
					f_r = (4*mass*np.pi*G*int_lm2[count]*np.real((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm1[count]*np.real(1j*m*(1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)

					#compute the phi force (theta in physics convention)
					phi_comb = (4*mass*np.pi*G*int_lm1[count]*np.real((1/(np.sqrt(2)))*(sph_harm_deriv(m,l,theta,phi)+(-1)**m*sph_harm_deriv(-m,l,theta,phi))))/(2*l+1)



					
					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)

			
				count += 1

			else:

				pass

	# sum each contribtuion term by term to get total forces
	F_r = np.sum(F_r, axis=0)
	F_theta = np.sum(F_theta, axis=0)
	F_phi = np.sum(F_phi, axis=0)


	#interpolate forces
	interpolate_F_r = interpolate.interp1d(bincenters[1:],F_r, kind='cubic')
	interpolate_F_theta = interpolate.interp1d(bincenters[1:],F_theta, kind='cubic')
	interpolate_F_phi = interpolate.interp1d(bincenters[1:],F_phi, kind='cubic')	



	# transform forces to cartesian coordinnates (identity from cover of griffiths (let phi --> theta and theta --> phi))

	F_x = (np.sin(phi)*np.cos(theta)*interpolate_F_r(r_position)) + (np.cos(theta)*np.cos(phi)*interpolate_F_phi(r_position)) + (-np.sin(theta)*interpolate_F_theta(r_position))
	
	F_y = (np.sin(phi)*np.sin(theta)*interpolate_F_r(r_position)) + (np.cos(phi)*np.sin(theta)*interpolate_F_phi(r_position)) + (np.cos(theta)*interpolate_F_theta(r_position))
	
	F_z = (np.cos(phi)*interpolate_F_r(r_position)) + (-np.sin(phi)*interpolate_F_phi(r_position))

	#returns fuctions of r for forces in spherical and only values at the particles position for cartesian forces

	return(F_x,F_y,F_z,F_r,F_theta,F_phi,)
	






# calculate potnetial and forces at same time (r_position, theta, phi, mass) correspond to the particle of interest
# int_lm1 is for potential and theta/phi forces ... int_lm2 is for radial force only

def MEX_potential_force(l_max, r_position, theta, phi, mass, bincenters , int_lm1, int_lm2):

	# make function spherical harmonic derivative with respect to phi  (theta physics convention)
	def sph_harm_deriv(m,l,theta,phi):
		d_phi = m*(1/np.tan(phi))*special.sph_harm(m,l,theta,phi) + np.sqrt((l-m)*(l+m+1))*np.e**(-1j*theta)*special.sph_harm(m+1,l,theta,phi)

		return(d_phi)
	

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

	F_r = []
	F_theta = []
	F_phi = []


	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:


				if m < 0 :
					
					# compute the potential
					#compute the potential term for this l and m comb
					phi_comb = (-4*np.pi*G*int_lm1[count]*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)
					# append all the contributions of phi into another list of the phi
					phi_contribution.append(phi_comb)
					

					#compute the radial force
					f_r = (mass*4*np.pi*G*int_lm2[count]*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm1[count]*np.real(1j*m*(1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)

					#compute the phi force (theta in physics convention)
					f_phi = (4*mass*np.pi*G*int_lm1[count]*np.real((1/(1j*np.sqrt(2)))*(sph_harm_deriv(-m,l,theta,phi)-(-1)**m*sph_harm_deriv(m,l,theta,phi))))/(2*l+1)

					

					# append all the contributions of forces into another list of the forces
					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)



				elif m == 0 :

					# compute the potential
					#compute the potential term for this l and m comb
					phi_comb =(-4*np.pi*G*int_lm1[count]*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)
					# append all the contributions of phi into another list of the phi
					phi_contribution.append(phi_comb)
					

					#compute the radial force
					f_r =(4*mass*np.pi*G*int_lm2[count]*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm1[count]*np.real(1j*m*special.sph_harm(m,l,theta,phi)))/(2*l+1)

					#compute the phi force (theta in physics convention)
					f_phi = (4*mass*np.pi*G*int_lm1[count]*np.real(sph_harm_deriv(m,l,theta,phi)))/(2*l+1)


					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)



				elif m > 0:
					
					# compute the potential

					#compute the potential term for this l and m comb
					phi_comb = (-4*np.pi*G*int_lm1[count]*np.real((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)
					# append all the contributions of phi into another list of the phi
					phi_contribution.append(phi_comb)
					

					#compute the radial force
					f_r = (4*mass*np.pi*G*int_lm2[count]*np.real((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm1[count]*np.real(1j*m*(1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)

					#compute the phi force (theta in physics convention)
					phi_comb = (4*mass*np.pi*G*int_lm1[count]*np.real((1/(np.sqrt(2)))*(sph_harm_deriv(m,l,theta,phi)+(-1)**m*sph_harm_deriv(-m,l,theta,phi))))/(2*l+1)



					
					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)

			
				count += 1

			else:

				pass

	# sum each contribtuion term by term to get total potential
	phi_contribution_sum = np.sum(phi_contribution, axis=0)
	## interpolate data to return potential at point wanted
	interpolate_potential = interpolate.interp1d(bincenters[1:],phi_contribution_sum, kind='cubic')
	
	# potential at particles position 
	PHI = interpolate_potential(r_position)


	# sum each contribtuion term by term to get total forces
	F_r = np.sum(F_r, axis=0)
	F_theta = np.sum(F_theta, axis=0)
	F_phi = np.sum(F_phi, axis=0)

	#interpolate forces
	interpolate_F_r = interpolate.interp1d(bincenters[1:],F_r, kind='cubic')
	interpolate_F_theta = interpolate.interp1d(bincenters[1:],F_theta, kind='cubic')
	interpolate_F_phi = interpolate.interp1d(bincenters[1:],F_phi, kind='cubic')	



	# transform forces to cartesian coordinnates (identity from cover of griffiths (let phi --> theta and theta --> phi))

	F_x = (np.sin(phi)*np.cos(theta)*interpolate_F_r(r_position)) + (np.cos(theta)*np.cos(phi)*interpolate_F_phi(r_position)) + (-np.sin(theta)*interpolate_F_theta(r_position))
	
	F_y = (np.sin(phi)*np.sin(theta)*interpolate_F_r(r_position)) + (np.cos(phi)*np.sin(theta)*interpolate_F_phi(r_position)) + (np.cos(theta)*interpolate_F_theta(r_position))
	
	F_z = (np.cos(phi)*interpolate_F_r(r_position)) + (-np.sin(phi)*interpolate_F_phi(r_position))

	#returns fuctions of r for forces in spherical and only values at the particles position for cartesian forces

	return(F_x,F_y,F_z,F_r,F_theta,F_phi,PHI,phi_contribution_sum)

	

	
	





