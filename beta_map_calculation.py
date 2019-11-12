'''
Created on 5 Sep. 2017

@author: jakubb
'''

import ROOT
import numpy as np
import sys
import os
import glob

# CT voxel size in mm
ct_vox_x = 0.683594 
ct_vox_y = 0.683594 
ct_vox_z = 1.2 

# Matrix dimension
dim_x = 512 
dim_y = 415 
dim_z = 252

# Isocenter offset (needed to calculate the matrix size)
off_x = -173.942 - 0.5*ct_vox_x
off_y = -104.871 - 0.5*ct_vox_y
off_z = -202.8 - 0.5*ct_vox_z


# ########################
# MAIN
# ########################


if __name__ == '__main__':

	'''

	Script for the automatic beta plus map calculation from patient irradiation
	based on Prometheus simulations

	'''	

	beta_map = np.zeros((dim_x, dim_y, dim_z))

	# Directory where all root files are stored
	results = '/home/baran/git/Simulations_GATE/patient_positrons/code_tests'
	
	os.chdir(results)

	ns2s = 1E9
	offbeam_s = 1200
	offbeam_ns = offbeam_s*ns2s
	annihil = 'annihil'
	ca = 0
	cna = 0
	root_files_number = 5000

	x_lim = dim_x * ct_vox_x
	y_lim = dim_y * ct_vox_y
	z_lim = dim_z * ct_vox_z

	# Loop entrance
	for i in range(1,root_files_number + 1):

		if os.path.exists('ps_patient'+str(i)+'.root'):			

			f = ROOT.TFile.Open('ps_patient'+str(i)+'.root')
			print "Processed map: {0}".format(i)
 		
			tree = f.Get("PhaseSpace")
			entries = tree.GetEntries()
			
			for i, entry in enumerate(tree):
				'''
				print annihil[0]
				print entry.CreatorProcess[0]
				print type(annihil)
				print type(entry.CreatorProcess)
				print id(annihil)
				print id(entry.CreatorProcess)
				print entry.CreatorProcess[0] == annihil[0]
				print '\n'
				'''
				if (entry.CreatorProcess[0] == annihil[0] and entry.Time > offbeam_ns):
					ca = ca+1

					# Calculating matrix indices corresponding to the production point
					temp_x = int((entry.X-off_x)/ct_vox_x)
					temp_y = int((entry.Y-off_y)/ct_vox_y)
					temp_z = int((entry.Z-off_z)/ct_vox_z)
	
					# Checking if the production point is at the edge of the CT FOV
					if (entry.X >= x_lim):
						temp_x -= 1
					if (entry.Y >= y_lim):
						temp_y -= 1
					if (entry.Z >= z_lim):
						temp_z -= 1
	
					beta_map[temp_x, temp_y, temp_z] += 1
			f.Close()

		else:
			print '{} root file is missing!'.format (i)


	print 'Number of all annihilations: {}'.format (int(ca/2))
	# Divide by 2 due to the two gammas are produced 
	beta_map = int(np.divide(beta_map,2))
	print 'Maximum number of annihilations in primary voxel: {}'.format (np.max(beta_map))

	mm = np.memmap('beta_map_final_20.raw', dtype='float32', mode='w+', order = 'F', shape=np.shape(beta_map))
	mm[:] = beta_map[:]


