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
import glob
import scipy.ndimage.filters as img_filters


back_to_back_number = 1000000000
number_format = 'float32'
bytes_per_pixel = '4'
scanner = 'dualhead_3x4'

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 400. 

# FOV size in mm (needed to calculate the matrix size)
phantom_x = 250. 
phantom_y = 250. 
phantom_z = 300. 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 2.5 
voxel_y = 2.5
voxel_z = 2.5 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)


phantom_dimx = int(phantom_x/voxel_x)
phantom_dimy = int(phantom_y/voxel_y)
phantom_dimz = int(phantom_z/voxel_z)

# ########################
# CALCULATE TRUE RECONSTRUCTED IMAGE
# ########################

def calculate_ideal_phantom (mean_value):

	phantom = np.zeros((mat_x, mat_y, mat_z))

	phantom_counts = mean_value# back_to_back_number / (phantom_dimx * phantom_dimy * phantom_dimz)

	ideal_phantom = np.full((phantom_dimx,phantom_dimy,phantom_dimz), phantom_counts)
	
	phantom [(mat_x/2-phantom_dimx/2) : (mat_x/2+phantom_dimx/2) ,  (mat_y/2-phantom_dimy/2) : (mat_y/2+phantom_dimy/2), (mat_z/2-phantom_dimz/2) : (mat_z/2+phantom_dimz/2)] = ideal_phantom
		
	return phantom.flatten(), ideal_phantom.flatten()




# ########################
# MAIN
# ########################


if __name__ == '__main__':


	'''

	Script for calculating the NRMSD to find the optimal number of iterations for the PET reconstruction
	
	'''	

#	phantom, ideal_phantom = calculate_ideal_phantom(back_to_back_number / (phantom_dimx * phantom_dimy * phantom_dimz))	
#	ideal_phantom = ideal_phantom


	activity_image_path = '/home/baran/git/Simulations_GATE/patient_positrons/beta_map_final_20_0ff_20_on_pet_grid.raw'

	phantom = np.fromfile(activity_image_path, dtype=number_format).reshape((160,160,200), order = 'F')

	if (scanner[-1] == '4'):
		phantom = phantom[:,:,20:180]
		
	phantom = phantom

	pet_images_path = '/home/baran/Desktop/castor_recons/'+scanner+'_patient/recon_it'
	pet_images_ext = '.img'
	
	iter_no = 20

	calibration_factor = 1

	iterr = np.linspace(1, iter_no, iter_no)
	nrmsd = np.zeros(iter_no)

	for i in range (1, iter_no+1):

		image_path = pet_images_path + str(i) + pet_images_ext
		temp_image = np.fromfile(image_path, dtype=number_format).reshape((mat_x,mat_y,mat_z), order = 'F')

		# ONLY FOR PHANTOM CONSIDERATIONS
		#temp_image = temp_image[(mat_x/2-phantom_dimx/2) : (mat_x/2+phantom_dimx/2) ,  (mat_y/2-phantom_dimy/2) : (mat_y/2+phantom_dimy/2), (mat_z/2-phantom_dimz/2) : (mat_z/2+phantom_dimz/2)]
		#temp_image = img_filters.gaussian_filter(temp_image,sigma=1).flatten()

		# FOR FULL IMAGE CONSIDERATIONS
#		temp_image = img_filters.gaussian_filter(temp_image,sigma=1).flatten()

		temp_image = temp_image*calibration_factor
#		if i == 1:
#			temp_value = np.sum(temp_image)/(temp_image.shape[0])
		#	print temp_value
		#	ideal_phantom = calculate_ideal_phantom(temp_value).flatten()	

#		temp_image = np.fromfile(image_path, dtype=number_format).reshape((mat_x,mat_y,mat_z), order = 'F')
#		temp_image = img_filters.gaussian_filter(temp_image,(2,2,1)).flatten()


		# ONLY FOR PHANTOM CONSIDERATIONS
#		diff = np.subtract(temp_image, ideal_phantom)
#		RMS = np.sqrt(np.sum(np.power(ideal_phantom,2)))

		# FOR FULL IMAGE CONSIDERATIONS
		diff = np.subtract(temp_image, phantom)
		RMS = np.sqrt(np.sum(np.power(phantom,2)))


		diff_power = np.sum(np.power(diff,2))
		RMSD = np.sqrt(diff_power)

		nrmsd[i-1] = RMSD/RMS
		print 'it: {};  nrmsd = {}'.format(i, RMSD/RMS)
		
		if i == 3:
			fig, axs = plt.subplots(2,2)

			im1 = axs[0,0].imshow(temp_image[:,:,70])
			axs[0,0].set_title("temp_image_70")
			fig.colorbar(im1, axs[0,0])
		
			im2 = axs[0,1].imshow(phantom[:,:,70])
			axs[0,1].set_title("phantom_70")
			fig.colorbar(im2, axs[0,1])
		
			im3 = axs[1,0].imshow(temp_image[:,:,74])
			axs[1,0].set_title("temp_image_74")
			fig.colorbar(im3, axs[1,0])
		
			im4 = axs[1,1].imshow(phantom[:,:,74])
			axs[1,1].set_title("phantom_74")
			fig.colorbar(im4, axs[1,1])
		
			fig.subplots_adjust(right=0.8)
			cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
			fig.colorbar(im1, cax=cbar_ax)
	
			plt.show()
	
		#print i


	#print nrmsd


	print 'Optimal number of iteration (min NRMSD): {0}, NRMSD: {1}'.format(np.argmin(nrmsd)+1, nrmsd[np.argmin(nrmsd)])

	plt.plot (iterr, nrmsd, 'ro')
#	plt.axis([0, iter_no, 0.5, 0.7])

	plt.title('NRMSD calculated based on cubic water phantom simulation with 1G primary back-to-back source uniformly distributed within the phantom \n(Gaitanis, Anastasios, et al. "PET image reconstruction: A stopping rule for the MLEM algorithm based on properties of the updating coefficients."\n Computerized Medical Imaging and Graphics 34.2 (2010): 131-141.)\nSETUP: '+scanner , fontweight="bold")

	plt.xlabel('Iteration number' , fontweight="bold")
	plt.ylabel('NRMSD' , fontweight="bold")

	plt.grid(True)

	#plt.show()



