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
import pandas as pd

number_of_beams = 50

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
fov_x = 400. 
fov_y = 400. 
#fov_z = 500. 
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


def scanner_name (sc):

	sc_capital = 'BARREL'

	if sc[0] == 'b' and sc[-1] == '2':
		sc_capital = 'BARREL2'

	if sc[0] == 'b' and sc[-1] == '3':
		sc_capital = 'BARREL3'

	if sc[0] == 'd' and sc[-3] == '1':
		sc_capital = 'DUALHEAD_1x6'

	if sc[0] == 'd' and sc[-3] == '2':
		sc_capital = 'DUALHEAD_2x6'

	if sc[0] == 'd' and sc[-1] == '4':
		sc_capital = 'DUALHEAD_3x4'
		fov_z = 400. 
	return sc_capital


def calculate_reconstructed_activity (image_path, constrained_voxels, sc):

	'''

	Method to process the reconstructed image to get the activity profile only within the phantom 
	or within the phantom and 1 cm before (at the entrance of the beam) if the extend_data flag is True

	3D Gaussian filtering is applied prior to the activity profile calculation as the reconstructed images are not smoothed  
	
	'''	
	mat_zz = 200
	fov_zz = 500. 
	if sc[0] == 'd' and sc[-1] == '4':
		mat_zz = 160
		fov_zz = 400. 


	reconstructed_activity = np.fromfile(image_path, dtype=number_format).reshape((mat_x,mat_y,mat_zz), order = 'F')

	gaussian_smoothing_factor = 3
	offset_voxels = 4
	reconstructed_activity = img_filters.gaussian_filter(reconstructed_activity,(gaussian_smoothing_factor))

	a=int((fov_x/2-phantom_x/2.)/voxel_x)+constrained_voxels
	b=int((fov_x/2+phantom_x/2.)/voxel_x)-constrained_voxels
	c=int((fov_y/2-phantom_y/2.)/voxel_y)+constrained_voxels
	d=int((fov_y/2+phantom_y/2.)/voxel_y)-constrained_voxels
	e=int((fov_zz/2-phantom_z/2.)/voxel_z)
	f=int((fov_zz/2+phantom_z/2.)/voxel_z)

	## EXTRACTING PHANTOM + 2 VOXELS FOR BEAM ENTRANCE FIT	

	phantom_reconstructed_activity = reconstructed_activity[ a:b, c:d, e-offset_voxels:f ]
	z_linspace = np.linspace(-phantom_z/2 + voxel_z/2 - offset_voxels*voxel_z, phantom_z/2-voxel_z/2, int(phantom_z/voxel_z)+offset_voxels)
	z_reconstructed_activity = phantom_reconstructed_activity.sum(axis=0).sum(axis=0)
		
	return z_reconstructed_activity, z_linspace, phantom_reconstructed_activity



if __name__ == '__main__':

	colors_list = ['darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver', 'darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver']
	line_type_list = ['.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v', '.', '>', 'P', 'x', '*', 'o', 'v']

	scanner = sys.argv[1]
	constrained_voxels = 0

	if scanner[0] == 'd' and scanner[-1] == '4':
		fov_z = 400.
		mat_z = 160

	sum_50_3D = np.zeros((mat_x,mat_y,mat_z))
	sum_50_z = np.zeros((mat_x,mat_y,mat_z))

	full_51_3D = np.zeros((mat_x,mat_y,mat_z))
	full_51_z = np.zeros((mat_x,mat_y,mat_z))

	linspace_z = np.zeros((1))



	print scanner

	# Output file name with all fit parameters stored - depend from the choosen option
	f_n = 'multi_energy_dose_' + scanner + '.txt'

#	'''
	
	with open(os.path.join(recon_dir, f_n), 'w') as sum_file:

		sum_file.write('#setup\tenter_fit\tu_enter_fit\tdistal_fit\tu_distal_fit\tz_dist_fit\tu_z_dist_fit\tz_dose_fit\tu_z_dose_fit\tz_distal_dose_fit\tu_z_distal_dose_fit')

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

		for i in range (1, number_of_beams+2):

			print i
			fig, ax1 = plt.subplots()
			fig.set_size_inches(14, 6)

			#Setup the paths
			multi_recon_dir = os.path.join(recon_dir, 'multi_'+scanner)
			temp_recon_dir = os.path.join(multi_recon_dir,'beam_'+str(i))
			image_path = os.path.join(temp_recon_dir,'recon_it5.img')

			# Caluculate the activity 
			reconstructed, linspace_recon, phantom_reconstructed_activity = calculate_reconstructed_activity(image_path, constrained_voxels, scanner)
	
			if i==1:
				sum_50_3D = phantom_reconstructed_activity
				sum_50_z = reconstructed
				linspace_z = linspace_recon
			elif i==51:
				full_51_3D = phantom_reconstructed_activity
				full_51_z = reconstructed
			else:
				sum_50_3D = sum_50_3D + phantom_reconstructed_activity
				sum_50_z = sum_50_z + reconstructed	

			
			# Normalize the activity profile
			reconstructed = (reconstructed-min(reconstructed))/(max(reconstructed)-min(reconstructed))

			# Plot the recontructed profile
			ax1.plot(linspace_recon, reconstructed, color=colors_list[i], marker = line_type_list[i], linewidth = 0.5)

			local_maxs = argrelextrema(reconstructed[0:80], np.greater)[0]
			while(local_maxs[0]<5):
				local_maxs = local_maxs[1:]
			while(local_maxs[-1]>58):
				local_maxs = local_maxs[0:-1]
			if (local_maxs[-1] < 48):
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
			#ax2.set_ylabel('Dose [A.U.]', fontweight="bold", color='magenta')#, fontsize = 20)
			ax2.set_ylabel('True activity [A.U.]', fontweight="bold", color='magenta')#, fontsize = 20)
			ax2.tick_params('y', colors='magenta', width=2)#, labelsize=20) 

			## Plot the dose fit
			ax2.plot(fit_x_dose, sigmoid(fit_x_dose, *popt_dose), color='black', linewidth = 1)

			# Fitting distal fall-off based on last maximum 
			fit_x_distal = linspace_recon[local_maxs[-1]:local_maxs[-1]+16]
			fit_y_distal = reconstructed[local_maxs[-1]:local_maxs[-1]+16]
			## Fitting
			popt_distal, pcov_distal = curve_fit(sigmoid, fit_x_distal, fit_y_distal, method='dogbox', bounds=([0.5, 0., -1.],[1.2, 40., 1.]))
			## Fit parameters sigma
			perr_distal = np.sqrt(np.diag(pcov_distal))
			## Plot the distal fit
			ax1.plot(fit_x_distal, sigmoid(fit_x_distal, *popt_distal), color='black', linewidth = 1)

			l1 = 'Reconstructed activity profile'#: ' + temp_dir
			l2 = 'Reconstructed activity distal \nfall-off position: {0:.2f} ({1:.2f})'.format(popt_distal[1], perr_distal[1])

			# Calculating difference between distal and dose/activity fall-offs (dfference at 50% high fall-offs)
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
			sum_file.write('\n' + scanner + '_' + str(i) + '\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(popt_distal[1], perr_distal[1], popt_dose[1], perr_dose[1], z_diff_dose, z_diff_dose_error))

			# Save the figure
			os.chdir(multi_recon_dir)

			distal_entrance_fit_file = 'beam_' + str(i) + '_activity_distal_entrance_fits'
			distal_fit_file = 'beam_' + str(i) + '_activity_distal_fits'
			
			plt.savefig(distal_fit_file+'.pdf', dpi = 600, bbox_inches='tight')
			
			plt.close()
			
	
#	'''
	
	
	# ANALYSIS
	fits_results = np.genfromtxt(os.path.join(recon_dir, f_n), delimiter='\t')

	os.chdir(recon_dir)

	data_hist = fits_results[0:-1, 5] # ACTIVITY DIFF
	data_mean = np.mean (data_hist)
	data_std = np.std (data_hist)
	
	sc = scanner_name(scanner)

	#print fits_results
	print 'Scanner: ' + sc
	#print 'ROI: {0}x{0} mm^2'.format(5-(constrained_voxels*2*voxel_x/10), 5-(constrained_voxels*2*voxel_y/10))
	print 'Mean: {0}'.format(data_mean)
	print 'Std: {0}'.format(data_std)
	textstr = '\n'.join((
		r'$\mu$=%.2f mm' % (data_mean),
		r'$\sigma$=%.2f mm' % (data_std)))

	sc = scanner_name(scanner)
	fig, ax = plt.subplots()
	ax.hist(data_hist)
	#plot_title = 'Scanner: ' + sc + u'   ROI: {0}x{1} cm\u00b2'.format(5-(constrained_voxels*2*voxel_x/10), 5-(constrained_voxels*2*voxel_y/10))
	#ax.set_title(plot_title, fontweight = 'bold')
	ax.set_xlabel('Difference between reconstructed activity and true activity profile [mm]')
	ax.set_ylabel('Frequency')
	ax.text(0.05, 0.85, textstr,transform=ax.transAxes)
	ax.grid(linestyle = ':', linewidth = 0.75, color = 'grey', axis = 'y')
	pdf_name = 'hist_' + sc + '.pdf' 
	plt.savefig(pdf_name, dpi = 600, bbox_inches='tight')	
	#plt.show()
	
#	'''

