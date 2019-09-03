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
number_of_beams = 7
simulations_repo = '/home/baran/git/Simulations_GATE'
recon_dir = '/home/baran/Desktop/castor_recons'

cnao_fov_x = 99.2 
cnao_fov_y = 220.8 
cnao_fov_z = 249.6 

analysis_type = 'var_energy'
#analysis_type = 'one_energy'

root_results = 'proton_beam_results.root'

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
voxel_x = 5. 
voxel_y = 5. 
voxel_z = 5. 

# VOXEL size in mm (scalling factor: mm/pixel)
dose_voxel_x = .5 
dose_voxel_y = .5 
dose_voxel_z = .5

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)

activity_dimx = int(phantom_x/voxel_x)
activity_dimy = int(phantom_y/voxel_y)
activity_dimz = int(phantom_z/voxel_z)

dose_path_dir = '/home/baran/git/Simulations_GATE/proton_barrel/'
dose_name = 'dose_phantom' 

dose_dimx = int(phantom_x/dose_voxel_x)
dose_dimy = int(phantom_y/dose_voxel_y)
dose_dimz = int(phantom_z/dose_voxel_z)

def calculate_dose ():

	final_dose = np.zeros((dose_dimx, dose_dimy, dose_dimz))	
	for i in range (0, 20):
		image_path = dose_path_dir + dose_name + str(i+1) + '-Dose.raw'
		temp_image = np.fromfile(image_path, dtype=number_format).reshape((dose_dimx,dose_dimy,dose_dimz), order = 'F')
		#print temp_image.shape
		final_dose = final_dose + temp_image

	z_dose_high = final_dose.sum(axis=0).sum(axis=0)
	z_linspace_high = np.linspace(-phantom_z/2 + dose_voxel_z/2, phantom_z/2-dose_voxel_z/2, dose_dimz)


	z_dose_low = np.zeros((int(phantom_z/voxel_z)))
	for i in range(0, z_dose_high.shape[0]):
		ind = int(i/(voxel_z/dose_voxel_z))
		z_dose_low[ind] += z_dose_high[i]
	z_linspace_low = np.linspace(-phantom_z/2 + voxel_z/2, phantom_z/2-voxel_z/2, activity_dimz)


	return z_dose_low, z_dose_high, z_linspace_low, z_linspace_high, final_dose

def calculate_reconstructed_activity (image_path):

	reconstructed_activity = np.fromfile(image_path, dtype=number_format).reshape((mat_x,mat_y,mat_z), order = 'F')
	reconstructed_activity = img_filters.gaussian_filter(reconstructed_activity,(1))
	#print reconstructed_activity.shape
	a=int((fov_x/2-phantom_x/2.)/voxel_x)
	b=int((fov_x/2+phantom_x/2.)/voxel_x)
	c=int((fov_y/2-phantom_y/2.)/voxel_y)
	d=int((fov_y/2+phantom_y/2.)/voxel_y)
	e=int((fov_z/2-phantom_z/2.)/voxel_z)
	f=int((fov_z/2+phantom_z/2.)/voxel_z)
	## ONLY PHANTOM	
	#phantom_reconstructed_activity = reconstructed_activity[ a:b, c:d, e:f ]
	z_linspace = np.linspace(-phantom_z/2 + voxel_z/2, phantom_z/2-voxel_z/2, activity_dimz)
#	print z_linspace
	## PHANTOM + 3 VOXELS FOR BEAM ENTRANCE FIT	
	phantom_reconstructed_activity = reconstructed_activity[ a:b, c:d, e-5:f ]
	z_linspace = np.linspace(-phantom_z/2 + voxel_z/2 - 5*voxel_z, phantom_z/2-voxel_z/2, activity_dimz+5)
#	print z_linspace

	z_reconstructed_activity = phantom_reconstructed_activity.sum(axis=0).sum(axis=0)

		
	return z_reconstructed_activity, z_linspace, phantom_reconstructed_activity


def define_castor_scanner ():

	castor_scanner = 'BARREL1'
	vox_x = 5.0
	vox_y = 5.0
	vox_z = 5.0
	dim_x = 80 
	dim_y = 80
	dim_z = 100

	if scanner == 'cnao_lso':
		castor_scanner = 'CNAO_LSO'
		vox_x = 1.6
		vox_y = 1.6
		vox_z = 1.6
		dim_x = int(cnao_fov_x/dim_x) 
		dim_y = int(cnao_fov_y/dim_y)
		dim_z = int(cnao_fov_z/dim_z)
	
	if scanner == 'barrel2':
		castor_scanner = 'BARREL2'

	if scanner == 'barrel3':
		castor_scanner = 'BARREL3'

	if scanner == 'dualhead_1x6':
		castor_scanner = 'DUALHEAD_1x6'

	if scanner == 'dualhead_2x6':
		castor_scanner = 'DUALHEAD_2x6'

	if scanner == 'dualhead_3x4':
		castor_scanner = 'DUALHEAD_3x4'
		dim_z = 80
	

	return castor_scanner, vox_x, vox_y, vox_z, dim_x, dim_y, dim_z



def sigmoid(x, a, b, c):
     y = a / (1 + np.exp(-c*(x-b)))
     return y


if __name__ == '__main__':

		
	
	if analysis_type == 'one_energy':
		number_of_beams = 7
		if scanner != 'barrel':
			print 'One energy reconsrtruction could be only done for the barrel setup!'
			sys.exit(0)
	
	for i in range(0, number_of_beams):

		## Setting the paths
		temp_dir = 'proton'+str(i+1)+'_'+scanner
		if analysis_type == 'one_energy':
			temp_dir = 'proton_'+scanner+'_'+str(i+1)

		gate_dir = os.path.join(simulations_repo, temp_dir)
		castor_dir = os.path.join(recon_dir, temp_dir)
		castor_dir_0 = ''
		if i == 0:
			castor_dir_0 = castor_dir

		print 'GATE dir: ' + gate_dir			
		print 'CASTOR dir: ' + castor_dir + '\n'			
		os.chdir(gate_dir)
		
		## Merging root files and scp the data
		if not os.path.exists('proton_beam_results.root'):
			shell_command = 'hadd proton_beam_results.root proton_beam_results*.root'
			print 'Merging .root files into one...'
			print shell_command
			p(shell_command, bufsize=0, shell=True)
			print '... done!\n'

		shell_command = 'scp proton_beam_results.root ' + castor_dir
		print 'Copying  proton_beam_results.root file to CASTOR recon dir...'
		print shell_command
		p(shell_command, bufsize=0, shell=True)
		print '... done!\n'




		## Setting the scanner type car
		[scanner_type, vox_x, vox_y, vox_z, dim_x, dim_y, dim_z] = define_castor_scanner()

		os.chdir(castor_dir)

		## Converting results file to castor format		
		shell_command = 'castor-GATERootToCastor -ots -src -TOF_reso 500 -i proton_beam_results.root -s ' + scanner_type + ' -vb 0 -m ~/git/Simulations_GATE/system_matrix_' + scanner + '/' + scanner + '.mac' + ' -o results'
		print 'Converting proton_beam_results.root file to the castor format...'
		print shell_command
		p(shell_command, bufsize=0, shell=True)
		print '... done!\n'

		print glob.glob('ac_*.hdr')
		## Performing the atten_sens recon if needed		
		if (not os.path.exists('atten_sens_it1.hdr') or not os.path.exists('atten_sens_it1.img')):
			
			if i == 0:
				attenuation_map_hdr = glob.glob('ac_*.hdr')
				if len(attenuation_map_hdr) == 1:	
					shell_command = 'castor-recon -df /home/shared/Castor_Sensitivity_LM/noTOF/' + scanner + '_df.Cdh -dim ' + str(dim_x) + ',' + str(dim_y) + ',' + str(dim_z) + ' -vox ' + str(vox_x) + ',' + str(vox_y) + ',' + str(vox_z) + ' -fout atten_sens -th 0 -it 1:1 -opti SENS -img ' + attenuation_map_hdr + ' -vb 2'
					print 'Recon merged atten_sens file ...'
					p(shell_command, bufsize=0, shell=True)
					print shell_command
					print '... done!\n'
				elif len(attenuation_map_hdr) == 0:
					print 'There is no attenuation map to perform the sensitivity and attenuation corrections merging...\n '
					print attenuation_map_hdr
					sys.exit(0)
				else :
					print 'There is TOO MUCH attenuation map to perform the sensitivity and attenuation corrections merging...\n '
					print attenuation_map_hdr
					sys.exit(0)
									
			else:
					shell_command = 'scp ' + castor_dir_0 + '/atten_sens* ./'
					print 'Copying  atten_sens files to CASTOR recon dir...'
					print shell_command
					p(shell_command, bufsize=0, shell=True)
					print '... done!\n'
				
			shell_command = 'scp proton_beam_results.root ' + castor_dir
			print 'Copying  proton_beam_results.root file to CASTOR recon dir...'
			print shell_command
			p(shell_command, bufsize=0, shell=True)
			print '... done!\n'

		## Performing the reconstruction		
		shell_command = 'castor-recon -df ./results_df.Cdh -dim ' + str(dim_x) + ',' + str(dim_y) + ',' + str(dim_z) + ' -vox ' + str(vox_x) + ',' + str(vox_y) + ',' + str(vox_z) + ' -fout recon -th 0 -it 10:1 -sens ./atten_sens_it1.hdr -vb 2' 
		print 'Reconstruction ...'
		print shell_command
		p(shell_command, bufsize=0, shell=True)
		print '... done!\n\n'
	'''


	## ONE ENRGY ANALYSIS


	colors_list = ['darkmagenta', 'red', 'green', 'blue', 'gold', 'peru', 'silver']
#	line_type_list = ['-', ':', '-.', '--', '*', 'o', 'x']
	line_type_list = ['.', '>', 'P', 'x', '*', 'o', 'v']
	
	if analysis_type == 'one_energy': 

		print 'Calculating dose ...'
		dose_low, dose_high, linspace_dose_low, linspace_dose_high, final_dose = calculate_dose()

		dose = dose_high
		x_dose = linspace_dose_high
	
		low = 0	
		if low:
			dose = dose_low
			x_dose = linspace_dose_low


	

		for i in range (0, number_of_beams):

			print i
			fig, ax1 = plt.subplots()

			dose = (dose-min(dose))/(max(dose)-min(dose))

			temp_dir = 'proton_'+scanner+'_'+str(i+1)
			castor_dir = os.path.join(recon_dir, temp_dir)	
			image_path = os.path.join(castor_dir,'recon_it7.img')
			reconstructed, linspace_recon, phantom_reconstructed_activity = calculate_reconstructed_activity(image_path)

			reconstructed = (reconstructed-min(reconstructed))/(max(reconstructed)-min(reconstructed))

			ax1.plot(linspace_recon, reconstructed, color=colors_list[i], marker = line_type_list[i], linewidth = 0.5)


			# FIT DISTAL 22 - 30

			local_maxs = argrelextrema(reconstructed[0:40], np.greater)[0]
			print local_maxs
			fit_x_distal = linspace_recon[local_maxs[-1]:local_maxs[-1]+6]
			fit_y_distal = reconstructed[local_maxs[-1]:local_maxs[-1]+6]

			popt_distal, pcov_distal = curve_fit(sigmoid, fit_x_distal, fit_y_distal, method='dogbox', bounds=([0.5, 0., -1.],[1.2, 40., 1.]))
			print(popt_distal)
			perr_distal = np.sqrt(np.diag(pcov_distal))
#			print(pcov)
#			print perr
			ax1.plot(fit_x_distal, sigmoid(fit_x_distal, *popt_distal), color='black', linewidth = 1)



			# FIT ENTRANCE 0 - 10
			fit_x_entrance = linspace_recon[local_maxs[0]-5:local_maxs[0]+1]
			fit_y_entrance = reconstructed[local_maxs[0]-5:local_maxs[0]+1]

			popt_entrance, pcov_entrance = curve_fit(sigmoid, fit_x_entrance, fit_y_entrance, method='dogbox', bounds=([-1., -120., -1.],[1., -80., 1.]))
			print(popt_entrance)
			perr_entrance = np.sqrt(np.diag(pcov_entrance))
#			print(pcov)
#			print perr
			ax1.plot(fit_x_entrance, sigmoid(fit_x_entrance, *popt_entrance), color='black', linewidth = 1)
		
			ax1.set_xlabel('Position [mm]', fontweight="bold")#, fontsize = 20)	
			ax1.set_ylabel('Activity [A.U.]', fontweight="bold")#, fontsize = 20)
			ax1.grid(linestyle = ':', linewidth = 0.75, color = 'grey')
	 		ax1.tick_params('x', width=2)#, labelsize=20)   
			ax1.tick_params('y', width=2)#, labelsize=20)   
	
			ax2 = ax1.twinx()    
			ax2.plot(x_dose, dose, color = 'magenta', linestyle = '-', linewidth = 3.)
			ax2.set_ylabel('Dose [A.U.]', fontweight="bold", color='magenta')#, fontsize = 20)
			ax2.tick_params('y', colors='magenta', width=2)#, labelsize=20) 
			l1 = 'proton'+str(i)
			l2 = 'distal z = {0:.2f} ({1:.2f})'.format(popt_distal[1], perr_distal[1])
			l3 = 'entrance z = {0:.2f} ({1:.2f})'.format(popt_entrance[1], perr_entrance[1])
			z_diff = popt_distal[1]-popt_entrance[1]
			z_diff_error = np.sqrt(perr_distal[1]*perr_distal[1] + perr_entrance[1]*perr_entrance[1])
			l4 = 'z_diff = {0:.2f} ({1:.2f})'.format(z_diff, z_diff_error)
			ax1.plot([], [], ' ')
			ax1.legend([l1, l2, l3, l4], loc='lower center')#, fontsize = 25)
			ax2.legend(['Deposited dose'], loc='upper right')#, fontsize = 25)


			print l4
			if i == 6:
				plt.show()
	'''

