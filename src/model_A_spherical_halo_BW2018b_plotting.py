# G=1 

# this code will make plots of positions and velocities of the model_A_spherical_halo_BW2018b

# importing useful packages 
import mpl_toolkits.mplot3d as mpl
import numpy as np
import matplotlib.pyplot as plt

# importing data for model_A_spherical_halo_BW2018b

loadfilepath = "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/halos/model_A_spherical_halo_BW2018b.ascii"
m,x,y,z,vx,vy,vz = np.loadtxt(loadfilepath,skiprows=1,unpack= True) 

# m is mass in units of 2.302e9 solar masses
# x,y,z is position in kpc 
# vx,vy,vz is velocity in units of 100 km/s

# make plots of x-y,x-z,y-z planes, and x-vx, y-vy, z-vz

# x-y plane
plt.hist2d(x, y,bins=50,cmap='Blues')
plt.ylabel("y [kpc]")
plt.xlabel("x [kpc]")
plt.title("x-y plane of model_A_spherical_halo_BW2018b")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fullFileName= "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/model_A_spherical_halo_BW2018b_plots/x-y_plane.eps"
plt.savefig(fullFileName)
plt.clf() 

#x-z plane
plt.hist2d(x, z,bins=50,cmap='Blues')
plt.ylabel("z [kpc]")
plt.xlabel("x [kpc]")
plt.title("x-z plane of model_A_spherical_halo_BW2018b")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fullFileName= "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/model_A_spherical_halo_BW2018b_plots/x-z_plane.eps"
plt.savefig(fullFileName)
plt.clf() 

#y-z plane 
plt.hist2d(y, z,bins=50,cmap='Blues')
plt.ylabel("y [kpc]")
plt.xlabel("z [kpc]")
plt.title("y-z plane of model_A_spherical_halo_BW2018b")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fullFileName= "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/model_A_spherical_halo_BW2018b_plots/y-z_plane.eps"
plt.savefig(fullFileName)
plt.clf() 

# x-vx
plt.hist2d(x, vx,bins=40,cmap='Blues')
plt.ylabel("vx [100 km/s]")
plt.xlabel("x [kpc]")
plt.title("x-vx position velocity distribution")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fullFileName= "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/model_A_spherical_halo_BW2018b_plots/x-vx_pos_vel_dist.eps"
plt.savefig(fullFileName)
plt.clf() 

# y-vy
plt.hist2d(y, vy,bins=40,cmap='Blues')
plt.ylabel("vy [100 km/s]")
plt.xlabel("y [kpc]")
plt.title("x-vx position velocity distribution")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fullFileName= "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/model_A_spherical_halo_BW2018b_plots/y-vy_pos_vel_dist.eps"
plt.savefig(fullFileName)
plt.clf() 

# z-vz
plt.hist2d(z, vz,bins=40,cmap='Blues')
plt.ylabel("vz [100 km/s]")
plt.xlabel("z [kpc]")
plt.title("z-vz position velocity distribution")
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
fullFileName= "/Users/Scott/Desktop/GitHub/CosmologicalPotentialExpansion/plots/model_A_spherical_halo_BW2018b_plots/z-vz_pos_vel_dist.eps"
plt.savefig(fullFileName)
plt.clf() 

