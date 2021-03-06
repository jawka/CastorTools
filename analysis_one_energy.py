import ROOT
import sys
import os
from multiprocessing import Pool
from subprocess import check_call as p
import matplotlib.pyplot as plt
import scipy.ndimage.filters as img_filters
from sklearn import preprocessing
import numpy as np
import glob
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema


scanner = 'barrel'
number_of_beams = 24
recon_dir = '/home/baran/Desktop/castor_recons'
dose_dir = '/home/shared/Gate/dose_maps/proton'
activity_dir = '/home/shared/Gate/activity_maps/proton'
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["font.weight"] = "bold"

number_format = 'float32'
bytes_per_pixel = '4'

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 

# FOV size in mm (needed to calculate the matrix size)
phantom_x = 50. 
phantom_y = 50. 
phantom_z = 200. 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 2.5 
voxel_y = 2.5 
voxel_z = 2.5 

# VOXEL size in mm (scalling factor: mm/pixel)
dose_voxel_x = .5 
dose_voxel_y = .5 
dose_voxel_z = .5

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)

def sigmoid(x, a, b, c):

	y = a / (1 + np.exp(-c*(x-b)))
	return y


def calculate_reconstructed_activity (image_path):

	'''

	Method to process the reconstructed image to get the activity profile only within the phantom 
	or within the phantom and 1 cm before (at the entrance of the beam) if the extend_data flag is True

	3D Gaussian filtering is applied prior to the activity profile calculation as the reconstructed images are not smoothed  
	
	'''	

	reconstructed_activity = np.fromfile(image_path, dtype=number_format).reshape((mat_x,mat_y,mat_z), order = 'F')

	gaussian_smoothing_factor = 3
	offset_voxels = 4
	reconstructed_activity = img_filters.gaussian_filter(reconstructed_activity,(gaussian_smoothing_factor))

	a=int((fov_x/2-phantom_x/2.)/voxel_x)
	b=int((fov_x/2+phantom_x/2.)/voxel_x)
	c=int((fov_y/2-phantom_y/2.)/voxel_y)
	d=int((fov_y/2+phantom_y/2.)/voxel_y)
	e=int((fov_z/2-phantom_z/2.)/voxel_z)
	f=int((fov_z/2+phantom_z/2.)/voxel_z)

	## EXTRACTING PHANTOM + 2 VOXELS FOR BEAM ENTRANCE FIT	

	phantom_reconstructed_activity = reconstructed_activity[ a:b, c:d, e-offset_voxels:f ]
	z_linspace = np.linspace(-phantom_z/2 + voxel_z/2 - offset_voxels*voxel_z, phantom_z/2-voxel_z/2, int(phantom_z/voxel_z)+offset_voxels)

	z_reconstructed_activity = phantom_reconstructed_activity.sum(axis=0).sum(axis=0)
		
	return z_reconstructed_activity, z_linspace, phantom_reconstructed_activity



if __name__ == '__main__':

	colors_list = ['darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue']
	line_type_list = ['.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x']
	
	# Output file name with all fit parameters stored - depend from the choosen option
	f_n = 'one_energy.txt'

	with open(os.path.join (recon_dir, f_n), 'w') as sum_file:

		sum_file.write('#setup\tdistal_fit\tu_distal_fit\tz_dose_fit\tu_z_dose_fit\tz_distal_dose_fit\tu_z_distal_dose_fit')

		print 'Calculating dose/true activity ...'		
		# Loading the precalculated dose profile
		#dose = np.load(os.path.join(dose_dir, 'dose_1D.npy'))
		dose = np.load(os.path.join(activity_dir, 'true_activity_1D.npy'))
		x_dose = np.linspace(-phantom_z/2 + dose_voxel_z/2, phantom_z/2-dose_voxel_z/2, int(phantom_z/dose_voxel_z))
		# Normalize the dose profile
		dose = (dose-min(dose))/(max(dose)-min(dose))
		# Fitting dose fall-off based on maximum 
		fit_x_dose = x_dose[np.argmax(dose) : np.argmax(dose)+25]
		fit_y_dose = dose[np.argmax(dose) : np.argmax(dose)+25]
		## Fitting
		popt_dose, pcov_dose = curve_fit(sigmoid, fit_x_dose, fit_y_dose, method='dogbox', bounds=([-1., 15., -1.],[1, 50., 1.]))
		## Fit parameters sigma
		perr_dose = np.sqrt(np.diag(pcov_dose))
		print popt_dose	

		for i in range (0, number_of_beams):

			print i
			fig, ax1 = plt.subplots()
			fig.set_size_inches(14, 6)

			#Setup the paths
			temp_dir = 'proton_'+scanner+'_'+str(i+1)
			castor_dir = os.path.join(recon_dir, temp_dir)	
			image_path = os.path.join(castor_dir,'recon_it7.img')
			if voxel_x == 2.5:
				image_path = os.path.join(castor_dir,'recon_2_5_it7.img')

			# Caluculate the activity 
			reconstructed, linspace_recon, phantom_reconstructed_activity = calculate_reconstructed_activity(image_path)
			
			# Normalize the activity profile
			reconstructed = (reconstructed-min(reconstructed))/(max(reconstructed)-min(reconstructed))

			# Plot the recontructed profile
			ax1.plot(linspace_recon, reconstructed, color=colors_list[i], marker = line_type_list[i], linewidth = 0.5)

			local_maxs = argrelextrema(reconstructed[0:80], np.greater)[0]
			#print local_maxs
			while(local_maxs[0]<5):
				local_maxs = local_maxs[1:]
			while(local_maxs[-1]>58):
				local_maxs = local_maxs[0:-1]

			print local_maxs
			# Figure's handling
			ax1.set_xlabel('Position [mm]', fontweight="bold")#, fontsize = 20)	
			ax1.set_ylabel('Activity [A.U.]', fontweight="bold")#, fontsize = 20)
			ax1.grid(linestyle = ':', linewidth = 0.75, color = 'grey')
	 		ax1.tick_params('x', width=2)#, labelsize=20)   
			ax1.tick_params('y', width=2)#, labelsize=20)   
	
			ax2 = ax1.twinx()    
			ax2.plot(x_dose, dose, color = 'magenta', linestyle = '-', linewidth = 3.)
			#ax2.set_ylabel('Dose [A.U.]', fontweight="bold", color='magenta')#, fontsize = 20)
			ax2.set_ylabel('True activity [A.U.]', fontweight="bold", color='magenta')#, fontsize = 20)
			ax2.tick_params('y', colors='magenta', width=2)#, labelsize=20) 

			## Plot the dose fit
			ax2.plot(fit_x_dose, sigmoid(fit_x_dose, *popt_dose), color='black', linewidth = 1)

			# Fitting distal fall-off based on last maximum 
			fit_x_distal = linspace_recon[local_maxs[-1]:local_maxs[-1]+12]
			fit_y_distal = reconstructed[local_maxs[-1]:local_maxs[-1]+12]
			## Fitting
			popt_distal, pcov_distal = curve_fit(sigmoid, fit_x_distal, fit_y_distal, method='dogbox', bounds=([0.5, 0., -1.],[1.2, 40., 1.]))
			## Fit parameters sigma
			perr_distal = np.sqrt(np.diag(pcov_distal))
			## Plot the distal fit
			ax1.plot(fit_x_distal, sigmoid(fit_x_distal, *popt_distal), color='black', linewidth = 1)
	
			l1 = 'Reconstructed activity profile'#: ' + temp_dir
			l2 = 'Reconstructed activity distal \nfall-off position: {0:.2f} ({1:.2f})'.format(popt_distal[1], perr_distal[1])
			
			# Calculating difference between distal and dose fall-offs (dfference at 50% high fall-offs)
			z_diff_dose = popt_dose[1]-popt_distal[1]
			# Calculating sigma of the difference between distal and dose fall-offs
			z_diff_dose_error = np.sqrt(perr_distal[1]*perr_distal[1] + perr_dose[1]*perr_dose[1])

			# Figure's handling
			f_size = 12	
			plt.rcParams["axes.labelweight"] = "bold"
			ax1.legend([l1, l2], loc='upper left')#, fontsize = 25)
			ax2.plot([], [], ' ')
			l5 = 'Fall-offs difference: {0:.2f} ({1:.2f})'.format(z_diff_dose, z_diff_dose_error)
			ax2.legend(['Production activity profile', 'Production activity distal fall-off \nposition: {0:.2f} ({1:.2f})'.format(popt_dose[1], perr_dose[1]), l5], loc='upper right', fontsize = f_size)#, fontsize = 25)

			# Writing to summary file
			sum_file.write('\n' + scanner + '_' + str(1+i) + '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(popt_distal[1], perr_distal[1], popt_dose[1], perr_dose[1], z_diff_dose, z_diff_dose_error))

			# Save the figure
			os.chdir(castor_dir)

			distal_fit_file = 'distal_fit_'

			plt.savefig(distal_fit_file+temp_dir+'.pdf', dpi = 1200, bbox_inches='tight')

			#plt.show()
