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

number_of_beams = 24
recon_dir = '/home/baran/Desktop/castor_recons'
dose_dir = '/home/shared/Gate/dose_maps/proton'
activity_dir = '/home/shared/Gate/activity_maps/proton'
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["font.weight"] = "bold"

number_format = 'float32'
bytes_per_pixel = '4'

plot_all = 0
extend_data = 1

# FOV size in mm (needed to calculate the matrix size)
phantom_x = 50. 
phantom_y = 50. 
phantom_z = 200. 

# VOXEL size in mm (scalling factor: mm/pixel)
dose_voxel_x = .5 
dose_voxel_y = .5 
dose_voxel_z = .5

vox_file = 2.5

def define_castor_scanner (gantry):

	'''

	Method to set the scanner parameters
	
	'''	

	castor_scanner = 'BARREL1'
	vox_x = 2.5
	vox_y = 2.5
	vox_z = 2.5
	dim_x = 160 
	dim_y = 160
	dim_z = 200
	
	if gantry == 'dualhead_3x4':
		castor_scanner = 'DUALHEAD_3x4'
		dim_z = 160
	
	if gantry == 'barrel2':
		castor_scanner = 'BARREL2'

	if gantry == 'barrel3':
		castor_scanner = 'BARREL3'

	if gantry == 'dualhead_1x6':
		castor_scanner = 'DUALHEAD_1x6'

	if gantry == 'dualhead_2x6':
		castor_scanner = 'DUALHEAD_2x6'

	
	return castor_scanner, vox_x, vox_y, vox_z, dim_x, dim_y, dim_z



def sigmoid(x, a, b, c):
     y = a / (1 + np.exp(-c*(x-b)))
     return y


def calculate_reconstructed_activity (image_path, gantry):

	'''

	Method to process the reconstructed image to get the activity profile only within the phantom 
	or within the phantom and 1 cm before (at the entrance of the beam) if the extend_data flag is True

	3D Gaussian filtering is applied prior to the activity profile calculation as the reconstructed images are not smoothed  
	
	'''	

	castor_scanner, vox_x, vox_y, vox_z, dim_x, dim_y, dim_z = define_castor_scanner (gantry)
	fov_x = vox_x*dim_x
	fov_y = vox_y*dim_y
	fov_z = vox_z*dim_z

	reconstructed_activity = np.fromfile(image_path, dtype=number_format).reshape((dim_x,dim_y,dim_z), order = 'F')
	gaussian_smoothing_factor = 3
	offset_voxels = 4
	reconstructed_activity = img_filters.gaussian_filter(reconstructed_activity,(gaussian_smoothing_factor))

	a=int((fov_x/2-phantom_x/2.)/vox_x)
	b=int((fov_x/2+phantom_x/2.)/vox_x)
	c=int((fov_y/2-phantom_y/2.)/vox_y)
	d=int((fov_y/2+phantom_y/2.)/vox_y)
	e=int((fov_z/2-phantom_z/2.)/vox_z)
	f=int((fov_z/2+phantom_z/2.)/vox_z)
	 
	## EXTRACTING PHANTOM + 2 VOXELS FOR BEAM ENTRANCE FIT	
	phantom_reconstructed_activity = reconstructed_activity[ a:b, c:d, e-offset_voxels:f ]
	z_linspace = np.linspace(-phantom_z/2 + vox_z/2 - offset_voxels*vox_z, phantom_z/2-vox_z/2, int(phantom_z/vox_z)+offset_voxels)

	z_reconstructed_activity = phantom_reconstructed_activity.sum(axis=0).sum(axis=0)
		
	return z_reconstructed_activity, z_linspace, phantom_reconstructed_activity



if __name__ == '__main__':

	colors_list = ['darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue']
	line_type_list = ['.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x']

	scanner_list = ['barrel', 'barrel2', 'barrel3', 'dualhead_1x6', 'dualhead_2x6', 'dualhead_3x4']

	# Output file name with all fit parameters stored - depend from the choosen option
	f_n = 'var_energy.txt'

	with open(os.path.join (recon_dir, f_n), 'w') as sum_file:

		sum_file.write('#setup\tscanner\tdistal_fit\tu_distal_fit\tz_dose_fit\tu_z_dose_fit\tz_distal_dose_fit\tu_z_distal_dose_fit')
		print 'Calculating dose/true activity ...'

		for beam_no in range (1,8):
		
			# Loading the precalculated dose profile
			#dose = np.load(os.path.join(dose_dir+str(beam_no), 'dose_1D.npy'))
			dose = np.load(os.path.join((activity_dir+str(beam_no)), 'true_activity_1D.npy'))
			x_dose = np.linspace(-phantom_z/2 + dose_voxel_z/2, phantom_z/2-dose_voxel_z/2, int(phantom_z/dose_voxel_z))
			# Normalize the dose profile
			dose = (dose-min(dose))/(max(dose)-min(dose))
			# Fitting dose fall-off based on maximum 
			fit_x_dose = x_dose[np.argmax(dose) : np.argmax(dose)+20]
			fit_y_dose = dose[np.argmax(dose) : np.argmax(dose)+20]
			## Fitting
			popt_dose, pcov_dose = curve_fit(sigmoid, fit_x_dose, fit_y_dose, method='dogbox', bounds=([-1., 10., -1.],[1, 60., 1.]))
			## Fit parameters sigma
			perr_dose = np.sqrt(np.diag(pcov_dose))
			print popt_dose	

			for ii, gantry in enumerate (scanner_list):

				print gantry
				castor_scanner, vox_x, vox_y, vox_z, dim_x, dim_y, dim_z = define_castor_scanner (gantry)
				fig, ax1 = plt.subplots()
				fig.set_size_inches(14, 6)

				#Setup the paths
				temp_dir = 'proton'+str(beam_no)+'_'+gantry
				castor_dir = os.path.join(recon_dir, temp_dir)	
				image_path = os.path.join(castor_dir,'recon_2_5_it5.img')

				# Caluculate the activity 
				reconstructed, linspace_recon, phantom_reconstructed_activity = calculate_reconstructed_activity(image_path, gantry)
			
				# Normalize the activity profile
				reconstructed = (reconstructed-min(reconstructed))/(max(reconstructed)-min(reconstructed))

				# Plot the recontructed profile
				ax1.plot(linspace_recon, reconstructed, color=colors_list[beam_no], marker = line_type_list[beam_no], linewidth = 0.5)
	
				local_maxs = argrelextrema(reconstructed[0:80], np.greater)[0]
				while(local_maxs[0]<3):
					local_maxs = local_maxs[1:]
				while(local_maxs[-1]>58):
					local_maxs = local_maxs[0:-1]
				if (gantry == 'barrel' and beam_no == 1):
					local_maxs[-1] = 50
				print local_maxs



				# Figure's handling
				ax1.set_xlabel('Position [mm]', fontweight="bold")#, fontsize = 20)	
				ax1.set_ylabel('Activity [A.U.]', fontweight="bold")#, fontsize = 20)
				ax1.grid(linestyle = ':', linewidth = 0.75, color = 'grey')
		 		ax1.tick_params('x', width=2)#, labelsize=20)   
				ax1.tick_params('y', width=2)#, labelsize=20)   
	
				ax2 = ax1.twinx()    
				ax2.plot(x_dose, dose, color = 'magenta', linestyle = '-', linewidth = 3.)
				ax2.set_ylabel('True activity [A.U.]', fontweight="bold", color='magenta')#, fontsize = 20)
				ax2.tick_params('y', colors='magenta', width=2)#, labelsize=20) 

				## Plot the dose fit
				ax2.plot(fit_x_dose, sigmoid(fit_x_dose, *popt_dose), color='black', linewidth = 1)
		
				# Fitting distal fall-off based on last maximum 
				fit_x_distal = linspace_recon[local_maxs[-1]:local_maxs[-1]+16]
				fit_y_distal = reconstructed[local_maxs[-1]:local_maxs[-1]+16]

				## Fitting
				popt_distal, pcov_distal = curve_fit(sigmoid, fit_x_distal, fit_y_distal, method='dogbox', bounds=([0.5, 20., -1.],[1.2, 40., 1.]))
				## Fit parameters sigma
				perr_distal = np.sqrt(np.diag(pcov_distal))
				## Plot the distal fit
				ax1.plot(fit_x_distal, sigmoid(fit_x_distal, *popt_distal), color='black', linewidth = 1)

				fit_x_entrance = linspace_recon[0:local_maxs[0]+1]
				fit_y_entrance = reconstructed[0:local_maxs[0]+1]

				l1 = 'Reconstructed activity profile'#: ' + temp_dir
				l2 = 'Reconstructed activity distal \nfall-off position: {0:.2f} ({1:.2f})'.format(popt_distal[1], perr_distal[1])
	
				# Calculating difference between distal and dose fall-offs (dfference at 50% high fall-offs)
				z_diff_dose = popt_dose[1]-popt_distal[1]
				# Calculating sigma of the difference between distal and dose fall-offs
				z_diff_dose_error = np.sqrt(perr_distal[1]*perr_distal[1] + perr_dose[1]*perr_dose[1])
	
				# Figure's handling
				f_size = 12	
				ax1.legend([l1, l2], loc='lower left', fontsize = f_size)#, fontsize = 25)
				ax2.plot([], [], ' ')
				l5 = 'Fall-offs difference: {0:.2f} ({1:.2f})'.format(z_diff_dose, z_diff_dose_error)
				ax2.legend(['Production activity profile', 'Production activity distal fall-off \nposition: {0:.2f} ({1:.2f})'.format(popt_dose[1], perr_dose[1]), l5], loc='upper right', fontsize = f_size)#, fontsize = 25)
	
				# Writing to summary file
				sum_file.write('\n' + 'proton'+str(beam_no)+'\t'+ gantry + '\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}'.format( popt_distal[1], perr_distal[1], \
						popt_dose[1], perr_dose[1], z_diff_dose, z_diff_dose_error))

				# Save the figure
				os.chdir(castor_dir)

				distal_fit_file = 'distal_fit_'
				plt.savefig(distal_fit_file+temp_dir+'.pdf', dpi = 1200, bbox_inches='tight')
	
				#plt.show()
