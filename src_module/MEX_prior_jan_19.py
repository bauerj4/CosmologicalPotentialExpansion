# this code calculates MEX potential
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special 
import scipy.integrate as integrate
import scipy.interpolate as interpolate


# NOTE
# F_x, F_y, F_z, Potential, are arrays over bincenters[1:]  *****

# function to load halo path and return values in spherical coordinates
# with theta and phi having correct phase 
 
# need to put in fudge factors on theta and phi positions for now *************

def LoadHalo(loadhalopath):

	mass,x,y,z,vx,vy,vz = np.loadtxt(loadhalopath,skiprows=1,unpack= True) 
	# need to calculate the radius each partcile is at 
	r = np.sqrt((x*x+y*y+z*z))  
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
def bin_r_log(radius,n,max_r):
	n1 = n + 1
	# max radial bin                                   
	max_bin = np.log10(max_r)
	# make logarithimically spaced radial bins
	linear_space = np.linspace(0,max_bin,n1)
	log_bins = 10**(linear_space)
	bincenters = 0.5*(log_bins[1:]+log_bins[:-1])
	# convert the r data to the corresponding bin center the ppoint falls in, only for r
	r_trunc = np.copy(radius)
	i = 0
	while i < len(bincenters):

		count = 0
		for x in radius: 

			if log_bins[i] <= x < log_bins[i+1]:

				r_trunc[count] = bincenters[i]


			else:
				pass 
		
			count += 1 		

		i += 1

	return(bincenters,r_trunc,log_bins)

# calculate the radial basis functions 
# uses true particle r not truncated r 
# listed in list as (0,0), (1,-1), (1,0),(1,1),(2,-2) .......
def basis_function_mod(l_max, bincenters, log_bins,r, r_trunc, mass, theta_pos, phi_pos):
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
	
						count=0
						for x in r_trunc:

							if x == bincenters[i]:

								term = (1/(r[count]**2*deltar[i]))*mass[count]*np.real((special.sph_harm(m,l,theta_pos[count],phi_pos[count])).conjugate())

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

								term = (1/(r[count]**2*deltar[i]))*mass[count]* np.real(((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta_pos[count],phi_pos[count])-(-1)**m*special.sph_harm(m,l,theta_pos[count],phi_pos[count]))).conjugate())
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

								term = (1/(r[count]**2*deltar[i]))*mass[count]* np.real(((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta_pos[count],phi_pos[count])+(-1)**m*special.sph_harm(-m,l,theta_pos[count],phi_pos[count]))).conjugate())

								density_cont.append(term)

								



							else:
								pass

							count += 1 

						density.append(sum(density_cont))

						
						i += 1 

					basis_functions.append(density)


	return(basis_functions)


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
# this int is for potential only
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

	# counter to get basis function for each l,m 
	count = 0

	# int_lm is the contributin for each l/m 
	int_lm = []

	#integrate functions term by term 
	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:
				
				
				# function for integral1 
				func1 = basis_functions[count]*bincenters**(l+2)
				# function for integral2 
				func2 = basis_functions[count]/(bincenters**(l-1))

				#counter is just for basis_func calling
				count +=1	


				int_1 = [0]
				int_2 = [0]

				for i in range(len(bincenters)-1):
					
					this_contrib_1 = integrate.simps([func1[i],func1[i+1]],[bincenters[i],bincenters[i+1]])
					this_contrib_2 = integrate.simps([func2[i],func2[i+1]],[bincenters[i],bincenters[i+1]])

					int_1 += [this_contrib_1 + int_1[-1]]
					int_2 += [this_contrib_2 + int_2[-1]]

				# need to do this to combine terms to find total int_lm
				int_1 = int_1[1:]
				int_2 = int_2[1:]
				bincenters1 = bincenters[1:]
				int_1 = np.array(int_1)
				int_2 = np.array(int_2)
				bincenters1 = np.array(bincenters1)
				
				# now modify int_2 to subtract integral from 0 --> r, making it r --> inf
				# append 0 at beggining of Array, and then remove last elemnt of int_2
				int_2_modified = np.insert(int_2,0,0) 
				int_2_modified = int_2_modified[:-1]
				
				# combine the terms now
				int_lm_cont = (1/(bincenters1**(l+1)))*int_1 + (bincenters1**(l))*(int_2[-1]-int_2_modified)

				int_lm.append(int_lm_cont)

			else:

				pass

			



	return(int_lm)


# this int_lm only works for the forces, int_lm1 radial, int_lm2 theta/phi
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

				#counter is just for basis_func calling
				count +=1	


				int_1 = [0]
				int_2 = [0]

				for i in range(len(bincenters)-1):
					
					this_contrib_1 = integrate.simps([func1[i],func1[i+1]],[bincenters[i],bincenters[i+1]])
					this_contrib_2 = integrate.simps([func2[i],func2[i+1]],[bincenters[i],bincenters[i+1]])

					int_1 += [this_contrib_1 + int_1[-1]]
					int_2 += [this_contrib_2 + int_2[-1]]

				# need to do this to combine terms to find total int_lm
				int_1 = int_1[1:]
				int_2 = int_2[1:]
				bincenters1 = bincenters[1:]
				int_1 = np.array(int_1)
				int_2 = np.array(int_2)
				bincenters1 = np.array(bincenters1)
				
				# now modify int_2 to subtract integral from 0 --> r, making it r --> inf
				# append 0 at beggining of Array, and then remove last elemnt of int_2
				int_2_modified = np.insert(int_2,0,0) 
				int_2_modified = int_2_modified[:-1]
				
				# combine the terms now
				# radial force
				int_lm_cont1 = (-(l+1)/(bincenters1**(l+2)))*int_1 + (l*bincenters1**(l-1))*(int_2[-1]-int_2_modified)
				# theta/phi forces
				int_lm_cont2 = (1/(bincenters1**(l+2)))*int_1 + (bincenters1**(l-1))*(int_2[-1]-int_2_modified)


				int_lm1.append(int_lm_cont1)
				int_lm2.append(int_lm_cont2)

			else:

				pass		

			

	return(int_lm1, int_lm2)


# returns inetgarls for both force and potentials , int_lm1 is for potential ... int_lm2 is for r force only .. int_lm3 if for theta/phi force
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
	int_lm3 = []

	# counter to get basis function for each l,m 
	count = 0

	for l in L:
		for m in M: 
		

			if np.absolute(m) <= l:
				
				
				# function for integral1 
				func1 = basis_functions[count]*bincenters**(l+2)
				# function for integral2 
				func2 = basis_functions[count]/(bincenters**(l-1))

				#counter is just for basis_func calling
				count +=1	


				int_1 = [0]
				int_2 = [0]

				for i in range(len(bincenters)-1):
					
					this_contrib_1 = integrate.simps([func1[i],func1[i+1]],[bincenters[i],bincenters[i+1]])
					this_contrib_2 = integrate.simps([func2[i],func2[i+1]],[bincenters[i],bincenters[i+1]])

					int_1 += [this_contrib_1 + int_1[-1]]
					int_2 += [this_contrib_2 + int_2[-1]]

				# need to do this to combine terms to find total int_lm
				int_1 = int_1[1:]
				int_2 = int_2[1:]
				bincenters1 = bincenters[1:]
				int_1 = np.array(int_1)
				int_2 = np.array(int_2)
				bincenters1 = np.array(bincenters1)
				
				# now modify int_2 to subtract integral from 0 --> r, making it r --> inf
				# append 0 at beggining of Array, and then remove last elemnt of int_2
				int_2_modified = np.insert(int_2,0,0) 
				int_2_modified = int_2_modified[:-1]
				
				
				# combine the terms now
				
				#force in r 
				int_lm_cont2 = (-(l+1)/(bincenters1**(l+2)))*int_1 + (l*bincenters1**(l-1))*(int_2[-1]-int_2_modified)
				# potential 
				int_lm_cont1 = (1/(bincenters1**(l+1)))*int_1 + (bincenters1**(l))*(int_2[-1]-int_2_modified)
				# force in theta/phi
				int_lm_cont3 = (1/(bincenters1**(l+2)))*int_1 + (bincenters1**(l-1))*(int_2[-1]-int_2_modified)



				int_lm2.append(int_lm_cont2)
				int_lm1.append(int_lm_cont1)
				int_lm3.append(int_lm_cont3)

			else:

				pass		

			
			

	return(int_lm1, int_lm2, int_lm3)





# calculate the potential now function 
# theta, phi: is where you want the potential gven for 
# l_max must be the same as basis_function and MEX_integral_lm 
#  (theta, phi) correspond to particle of interest 

def MEX_potential(l_max,theta, phi, bincenters, int_lm):
	
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



	return(phi_contribution_sum, phi_contribution)



# calculate the forces on the particles 
# int order is same as the function that calculates it 
# (theta, phi, mass) correspond to the particle of interest
# int_lm1 is for r force ... int_lm2 is for theta/phi force
def MEX_force(l_max, theta, phi, mass, bincenters , int_lm1, int_lm2):

	# make function spherical harmonic derivative with respect to phi  (theta physics convention)
	def sph_harm_deriv(m,l,theta,phi):
		
		if l == m : 
			
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


	# BT 2-122 1st ed. calculation (gradient of equation)

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
					f_r = (mass*4*np.pi*G*int_lm1[count]*np.real((1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm2[count]*np.real(1j*m*(1/(1j*np.sqrt(2)))*(special.sph_harm(-m,l,theta,phi)-(-1)**m*special.sph_harm(m,l,theta,phi))))/((2*l+1)*np.sin(phi))

					#compute the phi force (theta in physics convention)
					f_phi = (4*mass*np.pi*G*int_lm2[count]*np.real((1/(1j*np.sqrt(2)))*(sph_harm_deriv(-m,l,theta,phi)-(-1)**m*sph_harm_deriv(m,l,theta,phi))))/(2*l+1)

					

					# append all the contributions of forces into another list of the forces
					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)



				elif m == 0 :

					
					#compute the radial force
					f_r =(4*mass*np.pi*G*int_lm1[count]*np.real(special.sph_harm(m,l,theta,phi)))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm2[count]*np.real(1j*m*special.sph_harm(m,l,theta,phi)))/((2*l+1)*np.sin(phi))

					#compute the phi force (theta in physics convention)
					f_phi = (4*mass*np.pi*G*int_lm2[count]*np.real(sph_harm_deriv(m,l,theta,phi)))/(2*l+1)


					F_r.append(f_r)
					F_theta.append(f_theta)
					F_phi.append(f_phi)



				elif m > 0:
					
					
					#compute the radial force
					f_r = (4*mass*np.pi*G*int_lm1[count]*np.real((1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/(2*l+1)

					#compute the theta force (phi in physics convention)
					f_theta = (4*mass*np.pi*G*int_lm2[count]*np.real(1j*m*(1/(np.sqrt(2)))*(special.sph_harm(m,l,theta,phi)+(-1)**m*special.sph_harm(-m,l,theta,phi))))/((2*l+1)*np.sin(phi))

					#compute the phi force (theta in physics convention)
					phi_comb = (4*mass*np.pi*G*int_lm2[count]*np.real((1/(np.sqrt(2)))*(sph_harm_deriv(m,l,theta,phi)+(-1)**m*sph_harm_deriv(-m,l,theta,phi))))/(2*l+1)



					
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


	# transform forces to cartesian coordinnates (identity from cover of griffiths (let phi --> theta and theta --> phi))

	F_x = (np.sin(phi)*np.cos(theta)*F_r) + (np.cos(theta)*np.cos(phi)*F_phi )+ (-np.sin(theta)*F_theta)
	
	F_y = (np.sin(phi)*np.sin(theta)*F_r) + (np.cos(phi)*np.sin(theta)*F_phi) + (np.cos(theta)*F_theta)
	
	F_z = (np.cos(phi)*F_r) + (-np.sin(phi)*F_phi)

	#returns F_x, F_y, F_z over array Bincenetrs[1:]

	return(F_x,F_y,F_z)
	

