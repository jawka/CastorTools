'''
Created on 5 Sep. 2017

@author: jakubb
'''

import ROOT
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from multiprocessing import Pool
from matplotlib.colors import LogNorm
import glob


back_to_back_number = 100000000000
number_format = "float"
bytes_per_pixel = "4"
scanner = "barrel"
phantom_name = "pmma_5_5_20_cm_shifted"

sensitivity_binary_castor = "/home/baran/git/CastorTools/castor_sensitivity_maps/sensitivity_"+scanner+".raw"
sensitivity_header_castor = "/home/baran/git/CastorTools/castor_sensitivity_maps/sensitivity_"+scanner+".hdr"

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 5. 
voxel_y = 5. 
voxel_z = 5. 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)

# Annihilation number per voxel
prim_per_voxel = float(back_to_back_number)/(mat_x*mat_y*mat_z)

# Maps path
sens_maps_path= '/home/baran/git/CastorTools/castor_sensitivity_maps/barrel_data'



# ########################
# MAIN
# ########################


if __name__ == '__main__':




	print "FOV size: {0} x {1} x {2} mm3".format(fov_x, fov_y, fov_z)
	print "Voxel size: {0} x {1} x {2} (mm/pixel)3".format(voxel_x, voxel_y, voxel_z)
	print "Matrix size: {0} x {1} x {2}\n".format(mat_x, mat_y, mat_z)

	os.chdir(sens_maps_path)
	sensitivity_path = 'final_sensitivity_map.npy'	
	sensitivity = np.load(sensitivity_path)

	print "\nSensitivity matrix counts/voxel statistics:"	
	print "Sensitivity matrix min: {}".format(sensitivity.min())
	print "Sensitivity matrix max: {}\n".format(sensitivity.max())


	fig, axs = plt.subplots(2,3)

	im1 = axs[0,0].imshow(sensitivity[:,:,0], vmin=0.0001, vmax=0.04, cmap='jet')#, norm = LogNorm())
	axs[0,0].set_title("Transaxial (xy) plane  Z=-25. cm",fontweight="bold")
	axs[0,0].invert_yaxis()
	axs[0,0].set_xticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[0,0].set_yticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[0,0].set_xlabel("X position of the scanner [cm]")
	axs[0,0].set_ylabel("Y position of the scanner [cm]")
	#fig.colorbar(im1, axs[0,0])
	
	im2 = axs[0,1].imshow(sensitivity[:,:,25], vmin=0.0001, vmax=0.04, cmap='jet')#, norm = LogNorm())
	axs[0,1].set_title("Transaxial (xy) plane: Z=-12.5 cm",fontweight="bold")
	axs[0,1].invert_yaxis()
	axs[0,1].set_xticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[0,1].set_yticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[0,1].set_xlabel("X position of the scanner [cm]")
	axs[0,1].set_ylabel("Y position of the scanner [cm]")
	#fig.colorbar(im2, axs[0,1])

	im3 = axs[1,0].imshow(sensitivity[:,:,35], vmin=0.0001, vmax=0.04, cmap='jet')#, norm = LogNorm())
	axs[1,0].set_title("Transaxial (xy) plane: Z=-7.5 cm",fontweight="bold")
	axs[1,0].invert_yaxis()
	axs[1,0].set_xticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[1,0].set_yticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[1,0].set_xlabel("X position of the scanner [cm]")
	axs[1,0].set_ylabel("Y position of the scanner [cm]")
	#fig.colorbar(im3, axs[1,0])

	im4 = axs[1,1].imshow(sensitivity[:,:,50], vmin=0.0001, vmax=0.04, cmap='jet')#, norm = LogNorm())
	axs[1,1].set_title("Transaxial (xy) plane: Z=0. cm",fontweight="bold")
	axs[1,1].invert_yaxis()
	axs[1,1].set_xticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[1,1].set_yticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[1,1].set_xlabel("X position of the scanner [cm]")
	axs[1,1].set_ylabel("Y position of the scanner [cm]")
	#fig.colorbar(im4, axs[1,1])

	im5 = axs[0,2].imshow(sensitivity[39,:,:], vmin=0.0001, vmax=0.04, cmap='jet')#, norm = LogNorm())
	axs[0,2].set_title("Axial (zy) plane: X = 0.0 cm",fontweight="bold")
	axs[0,2].invert_yaxis()
	axs[0,2].set_xticklabels(["","-25","-15","-5","5","15"])
	axs[0,2].set_yticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[0,2].set_xlabel("Z position of the scanner [cm]")
	axs[0,2].set_ylabel("Y position of the scanner [cm]")
	#fig.colorbar(im5, axs[0,2])

	im6 = axs[1,2].imshow(sensitivity[:,39,:], vmin=0.0001, vmax=0.04, cmap='jet')#, norm = LogNorm())
	axs[1,2].set_title("Axial (zx) plane: Y=0.0 cm",fontweight="bold")
	axs[1,2].invert_yaxis()
	axs[1,2].set_xticklabels(["","-25","-15","-5","5","15"])
	axs[1,2].set_yticklabels(["","-20","-15","-10","-5","0","5","10","15","20"])
	axs[1,2].set_xlabel("Z position of the scanner [cm]")
	axs[1,2].set_ylabel("X position of the scanner [cm]")
	#fig.colorbar(im6, axs[1,2])

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	fig.colorbar(im4, cax=cbar_ax)

	plt.show()



#	save_sensitivity_to_binary_argument(sensitivity)
#	save_sensitivity_to_header()



