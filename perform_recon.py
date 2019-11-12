import ROOT
import sys
import os
from multiprocessing import Pool
from subprocess import check_call as p
import matplotlib.pyplot as plt
from sklearn import preprocessing
import numpy as np
import glob


scanner = 'barrel'
number_of_beams = 7
simulations_repo = '/home/baran/git/Simulations_GATE'
recon_dir = '/home/baran/Desktop/castor_recons'

cnao_fov_x = 99.2 
cnao_fov_y = 220.8 
cnao_fov_z = 249.6 

#analysis_type = 'var_energy'
analysis_type = 'one_energy'

root_results = 'proton_beam_results.root'

number_format = 'float32'
bytes_per_pixel = '4'

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 


def define_castor_scanner ():

	castor_scanner = 'BARREL1'
	'''
	vox_x = 2.5
	vox_y = 2.5
	vox_z = 2.5
	dim_x = 160
	dim_y = 160
	dim_z = 200
	'''
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
		#dim_z = 200
		dim_z = 100

	if scanner == 'dualhead_3x4':
		castor_scanner = 'DUALHEAD_3x4'
		#dim_z = 160
		dim_z = 80

	return castor_scanner, vox_x, vox_y, vox_z, dim_x, dim_y, dim_z



def sigmoid(x, a, b, c):
     y = a / (1 + np.exp(-c*(x-b)))
     return y


if __name__ == '__main__':

		
	
	if analysis_type == 'one_energy':
		number_of_beams = 24
		if scanner != 'barrel':
			print 'One energy reconsrtruction could be only done for the barrel setup!'
			sys.exit(0)

	castor_dir_0 = ''
	
	for i in range(0, number_of_beams):

		## Setting the paths
		temp_dir = 'proton'+str(i+1)+'_'+scanner
		if analysis_type == 'one_energy':
			temp_dir = 'proton_'+scanner+'_'+str(i+1)

		gate_dir = os.path.join(simulations_repo, temp_dir)
		castor_dir = os.path.join(recon_dir, temp_dir)
		if i == 0:
			castor_dir_0 = castor_dir

		print 'GATE dir: ' + gate_dir			
		print 'CASTOR dir: ' + castor_dir + '\n'			
		os.chdir(gate_dir)
		'''
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
		'''
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
		if (not os.path.exists('atten_sens_2_5_it1.hdr') or not os.path.exists('atten_sens_2_5_it1.img')):
			
			if i == 0:
				attenuation_map_hdr = glob.glob('ac_*.hdr')
				if len(attenuation_map_hdr) == 1:	
					shell_command = 'castor-recon -df /home/shared/Castor_Sensitivity_LM/noTOF/' + scanner + '_df.Cdh -dim ' + str(dim_x) + ',' + str(dim_y) + ',' + str(dim_z) + ' -vox ' + str(vox_x) + ',' + str(vox_y) + ',' + str(vox_z) + ' -fout atten_sens_2_5 -th 0 -it 1:1 -opti SENS -img ' + attenuation_map_hdr + ' -vb 2'
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
				
		## Performing the reconstruction		
		shell_command = 'castor-recon -df ./results_df.Cdh -dim ' + str(dim_x) + ',' + str(dim_y) + ',' + str(dim_z) + ' -vox ' + str(vox_x) + ',' + str(vox_y) + ',' + str(vox_z) + ' -fout recon -th 0 -it 10:1 -sens ./atten_sens_it1.hdr -vb 2' 
		print 'Reconstruction ...'
		print shell_command
		p(shell_command, bufsize=0, shell=True)
		print '... done!\n\n'



