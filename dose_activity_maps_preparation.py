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
from sklearn import preprocessing

scanner = 'barrel3'
number_format = 'float32'
bytes_per_pixel = '4'
dose_name = 'dose_phantom' 
simulations_repo = '/home/baran/git/Simulations_GATE'
activity_dir 	= '/home/shared/Gate/activity_maps'
dose_dir 	= '/home/shared/Gate/dose_maps'
number_of_beams = 8

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 

cnao_fov_x = 99.2 
cnao_fov_y = 220.8 
cnao_fov_z = 249.6 

# FOV size in mm (needed to calculate the matrix size)
phantom_x = 50. 
phantom_y = 50. 
phantom_z = 200. 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 5. 
voxel_y = 5. 
voxel_z = 5. 

cnao_voxel_x = 1.6 
cnao_voxel_y = 1.6 
cnao_voxel_z = 1.6 

# VOXEL size in mm (scalling factor: mm/pixel)
dose_voxel_x = .5 
dose_voxel_y = .5 
dose_voxel_z = .5

# VOXEL size in mm (scalling factor: mm/pixel)
cnao_dose_voxel_x = 1.6 
cnao_dose_voxel_y = 1.6 
cnao_dose_voxel_z = 1.6

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)

cnao_mat_x = int(cnao_fov_x/cnao_voxel_x) 
cnao_mat_y = int(cnao_fov_y/cnao_voxel_y)
cnao_mat_z = int(cnao_fov_z/cnao_voxel_z)

activity_dimx = int(phantom_x/voxel_x)
activity_dimy = int(phantom_y/voxel_y)
activity_dimz = int(phantom_z/voxel_z)

dose_dimx = int(phantom_x/dose_voxel_x)
dose_dimy = int(phantom_y/dose_voxel_y)
dose_dimz = int(phantom_z/dose_voxel_z)

# ########################
# CALCULATE THE DOSE
# ########################

def calculate_dose ():

	for i in range (0, number_of_beams):

		final_dose = np.zeros((dose_dimx, dose_dimy, dose_dimz))	

		print 'Case: {0}'.format(i)

		for j in range (0, 20):
			f_name = 'proton' + str(i) + '_'+scanner

			if  i==0:
				f_name = 'proton_'+scanner
			dose_f_name = dose_name + str(j+1) + '-Dose.raw'
			image_path = os.path.join(simulations_repo, f_name, dose_f_name)
			#print image_path	
			temp_image = np.fromfile(image_path, dtype=number_format).reshape((dose_dimx,dose_dimy,dose_dimz), order = 'F')
			#print temp_image.shape
			final_dose = final_dose + temp_image
	
		z_dose = final_dose.sum(axis=0).sum(axis=0)
#		print z_dose.shape
		ff_name = 'proton' + str(i)
		if  i==0:
			ff_name = 'proton'

		final_dose_path = os.path.join(dose_dir, ff_name, 'dose_3D')
		z_dose_path = os.path.join(dose_dir, ff_name, 'dose_1D')

		#print final_dose_path
		#print z_dose_path	
			
		final_dose.flatten('F').astype('float32').tofile(final_dose_path+'.bin')
		z_dose.flatten('F').astype('float32').tofile(z_dose_path+'.bin')
		np.save(final_dose_path, final_dose)
		np.save(z_dose_path, z_dose)


# ########################
# CALCULATE THE TRUE ACTIVITY
# ########################



def calculate_true_activity ():

	
	no_root_files = 20;



	for i in range (0, 8):
		
		print i
		beta_map = np.zeros((dose_dimx, dose_dimy, dose_dimz))
		activity_dose_path = os.path.join(simulations_repo, 'proton' + str(i) + '_barrel')
		if (i == 0):
			activity_dose_path = os.path.join(simulations_repo, 'proton_barrel')
		#print activity_dose_path
		os.chdir(activity_dose_path)

		maps_list = glob.glob('positron_high*-Stop.img')
		if (i == 0):
			maps_list = glob.glob('positron_map_high*.npy')

		for f in maps_list:
			ff = ''
			if i ==0:
				ff = np.load(f) 
			else:
				ff = np.fromfile(f, dtype=number_format).reshape((dose_dimx,dose_dimy,dose_dimz), order = 'F')
			beta_map = beta_map + ff

		beta_map_z = beta_map.sum(axis=0).sum(axis=0)		
		p = 'proton'
		positron_path_3D = os.path.join(activity_dir, p, 'true_activity_3D')
		positron_path_1D = os.path.join(activity_dir, p, 'true_activity_1D')
		positron_path_3D_bin = os.path.join(activity_dir, p, 'true_activity_3D.bin')
		positron_path_1D_bin = os.path.join(activity_dir, p, 'true_activity_1D.bin')

		if (i != 0):
			positron_path_3D = os.path.join(activity_dir, p+str(i), 'true_activity_3D')
			positron_path_1D = os.path.join(activity_dir, p+str(i), 'true_activity_1D')
			positron_path_3D_bin = os.path.join(activity_dir, p+str(i), 'true_activity_3D.bin')
			positron_path_1D_bin = os.path.join(activity_dir, p+str(i), 'true_activity_1D.bin')

#		print np.shape(beta_map)
#		print np.shape(beta_map_z)
#		print np.max(beta_map)
#		print np.max(beta_map_z)
#		print positron_path_3D
#		print positron_path_1D

		beta_map.flatten('F').astype('int16').tofile(positron_path_3D_bin)
		beta_map_z.flatten('F').astype('int16').tofile(positron_path_1D_bin)
		np.save(positron_path_3D, beta_map)
		np.save(positron_path_1D, beta_map_z)

		

# ########################
# MAIN
# ########################


if __name__ == '__main__':



	print 'Calculating dose ...'
	calculate_dose()
	print '.. done!'	

	print 'Calculating true beta plus profile ...'
	calculate_true_activity()
	print '.. done!'	

