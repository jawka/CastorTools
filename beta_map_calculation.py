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
#ct_vox_x = 0.683594 
#ct_vox_y = 0.683594 
#ct_vox_z = 1.2 

ct_vox_x = 2.5 
ct_vox_y = 2.5 
ct_vox_z = 2.5 

# Matrix dimension
#dim_x = 512 
#dim_y = 415 
#dim_z = 252

dim_x = 160 
dim_y = 160 
dim_z = 200

# Isocenter offset (needed to calculate the matrix size)
#off_x = -173.942 - 0.5*ct_vox_x
#off_y = -104.871 - 0.5*ct_vox_y
#off_z = -202.8 - 0.5*ct_vox_z

off_x = -ct_vox_x*dim_x/2.
off_y = -ct_vox_y*dim_y/2.
off_z = -ct_vox_z*dim_z/2.

# ########################
# MAIN
# ########################

def save_activity_map_to_header (f_name):

	with open(f_name + '.hdr', 'w+') as f:
		f.write("!name of data  file :=  {0}\n".format(f_name+".raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("image data  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(dim_x))
		f.write("!matrix  size  [2] := {0}\n".format(dim_y))
		f.write("!matrix  size  [3] := {0}\n".format(dim_z))
		f.write("!number  format  := float32\n")
		f.write("!number  of bytes  per  pixel  := 4\n")
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(ct_vox_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(ct_vox_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(ct_vox_z))
		f.write("image  duration (sec) := 1\n")



if __name__ == '__main__':

	'''

	Script for the automatic beta plus map calculation from patient irradiation
	based on Prometheus simulations

	'''	

	beta_map_all = np.zeros((dim_x, dim_y, dim_z))
	beta_map_20_off_20_on = np.zeros((dim_x, dim_y, dim_z))
	beta_map_15_off_20_on = np.zeros((dim_x, dim_y, dim_z))

	# Directory where all root files are stored
	results = '/home/baran/git/Simulations_GATE/patient_positrons/prometheus_data'
		
	os.chdir(results)

	ns2s = 1E9
	offbeam_20_s = 1200
	offbeam_15_s = 900
	offbeam_20_ns = offbeam_20_s*ns2s
	offbeam_15_ns = offbeam_15_s*ns2s
	annihil = 'annihil'
	c_all = 0
	c_20_20 = 0
	c_15_20 = 0
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


				if (entry.CreatorProcess[0] == annihil[0] and entry.CreatorProcess[1] == annihil[1]):
					c_all = c_all + 1
					temp_x = int((entry.X-off_x)/ct_vox_x)
					temp_y = int((entry.Y-off_y)/ct_vox_y)
					temp_z = int((entry.Z-off_z)/ct_vox_z)

					if (entry.X >= x_lim):
						temp_x -= 1
					if (entry.Y >= y_lim):
						temp_y -= 1
					if (entry.Z >= z_lim):
						temp_z -= 1
	
					beta_map_all[temp_x, temp_y, temp_z] += 1

					if (entry.Time > offbeam_20_ns and entry.Time < 2*offbeam_20_ns):
						beta_map_20_off_20_on[temp_x, temp_y, temp_z] += 1
						c_20_20 = c_20_20 + 1

					if (entry.Time > offbeam_15_ns and entry.Time < offbeam_20_ns+offbeam_15_ns):
						beta_map_15_off_20_on[temp_x, temp_y, temp_z] += 1
						c_15_20 = c_15_20 + 1
			f.Close()

		else:
			print '{} root file is missing!'.format (i)

	os.chdir('/home/baran/git/Simulations_GATE/patient_positrons/')

	
	print '\nCase: ALL'
	print 'Number of all annihilations: {}'.format (int(c_all/2))
	# Divide by 2 due to the two gammas are produced 
	beta_map_all = np.divide(beta_map_all,2)
	beta_map_all = beta_map_all.astype(int)
	print 'Maximum number of annihilations in primary voxel: {}'.format (np.max(beta_map_all))

	mm_all = np.memmap('beta_map_final_all_pet_grid.raw', dtype='float32', mode='w+', order = 'F', shape=np.shape(beta_map_all))
	mm_all[:] = beta_map_all[:]



        print '\nCase: 20_20'
        print 'Number of all annihilations: {}'.format (int(c_20_20/2))
        # Divide by 2 due to the two gammas are produced 
        beta_map_20_off_20_on = np.divide(beta_map_20_off_20_on,2)
        beta_map_20_off_20_on = beta_map_20_off_20_on.astype(int)
        print 'Maximum number of annihilations in primary voxel: {}'.format (np.max(beta_map_20_off_20_on))

        mm_20_20 = np.memmap('beta_map_final_20_0ff_20_on_pet_grid.raw', dtype='float32', mode='w+', order = 'F', shape=np.shape(beta_map_20_off_20_on))
        mm_20_20[:] = beta_map_20_off_20_on[:]


        print '\nCase: 15_20'
        print 'Number of all annihilations: {}'.format (int(c_15_20/2))
        # Divide by 2 due to the two gammas are produced 
        beta_map_15_off_20_on = np.divide(beta_map_15_off_20_on,2)
        beta_map_15_off_20_on = beta_map_15_off_20_on.astype(int)
        print 'Maximum number of annihilations in primary voxel: {}'.format (np.max(beta_map_15_off_20_on))

        mm_15_20 = np.memmap('beta_map_final_15_0ff_20_on_pet_grid.raw', dtype='float32', mode='w+', order = 'F', shape=np.shape(beta_map_15_off_20_on))
        mm_15_20[:] = beta_map_15_off_20_on[:]

	


	save_activity_map_to_header('beta_map_final_all_pet_grid')
	save_activity_map_to_header('beta_map_final_20_0ff_20_on_pet_grid')
	save_activity_map_to_header('beta_map_final_15_0ff_20_on_pet_grid')




