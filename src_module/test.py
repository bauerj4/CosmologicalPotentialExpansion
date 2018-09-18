import MEX as MEX

loadhalopath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/halos/model_A_spherical_halo_BW2018b.ascii"
mass, r, theta_pos, phi_pos = MEX.LoadHalo(loadhalopath)

bincenters, r_trunc, log_bin = MEX.bin_r_log(r,32,1000)

basis_func = MEX.basis_function(0, bincenters, log_bin, r_trunc, mass, theta_pos, phi_pos)

MEX.MEX_potential_plot(0, 0 , 0, bincenters, basis_func)