import MEX as MEX
import numpy as np
import pylab as pl

loadhalopath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/halos/model_A_spherical_halo_BW2018b.ascii"
mass, r, theta_pos, phi_pos = MEX.LoadHalo(loadhalopath)

bincenters, r_trunc, log_bin = MEX.bin_r_log(r,32,1000)

basis_func = MEX.basis_function(0, bincenters, log_bin, r_trunc, mass, theta_pos, phi_pos)

int_lm = MEX.MEX_integral_lm_potential(0,basis_func,bincenters)

int_lm11, int_lm2 = MEX.MEX_integral_lm_force_potential(0,basis_func,bincenters)
print()
print(int_lm)
print()
print(int_lm11)



PHI, phi_contribution_sum, phi_contribution = MEX.MEX_potential(0, 11, 0, 0, bincenters, int_lm11)
PHI, phi_contribution_sum1, phi_contribution = MEX.MEX_potential(0, 11, 0, 0, bincenters, int_lm)
pl.plot(bincenters[1:], phi_contribution_sum)
pl.plot(bincenters[1:], phi_contribution_sum1)
pl.show()