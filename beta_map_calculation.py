'''
Created on 5 Sep. 2017

@author: jakubb
'''

import ROOT
import numpy as np
import sys
import os
#import matplotlib.pyplot as plt
#from multiprocessing import Pool
import glob


# Isocenter offset (needed to calculate the matrix size)
off_x = -173.942 
off_y = -104.871
off_z = -202.8

# CT voxel size in mm
ct_vox_x = 0.683594 
ct_vox_y = 0.683594 
ct_vox_z = 1.2 

# Matrix dimension
dim_x = 512 
dim_y = 415 
dim_z = 252




# ########################
# MAIN
# ########################


if __name__ == '__main__':

	beta_map = np.zeros((dim_x, dim_y, dim_z))

	results = '/home/baran/git/Simulations_GATE/patient_positrons/code_tests'
	
	os.chdir(results)


	ns2s = 1E9
	offbeam_s = 1200
	offbeam_ns = offbeam_s*ns2s
	annihil = 'annihil'
	ca = 0
	cna = 0

	for i in range(1,21):
		if os.path.exists('ps_patient'+str(i)+'.root'):			
		#	try:
			f = ROOT.TFile.Open('ps_patient'+str(i)+'.root')
			print "Processed map: {0}".format(i)
		#	except:
		#		print "File {0} doesn't open properly!".format('ps_patient'+str(i)+'.root')
		#		sys.exit(2)
 		
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
					temp_x = int((entry.X-off_x)/ct_vox_x)
					temp_y = int((entry.Y-off_y)/ct_vox_y)
					temp_z = int((entry.Z-off_z)/ct_vox_z)
	
					if (entry.X >= 350.000128):
						temp_x -= 1
					if (entry.Y >= 283.69151):
						temp_y -= 1
					if (entry.Z >= 302.4):
						temp_z -= 1
	
					beta_map[temp_x, temp_y, temp_z] += 1
			f.Close()

		else:
			print '{} root file is missing!'.format (i)


	print 'Number of all annihilations: {}'.format (int(ca/2))
	beta_map = np.divide(beta_map,2)
	print 'Maximum number of annihilations in primary voxel: {}'.format (np.max(beta_map))

	mm = np.memmap('beta_map_final_20.raw', dtype='float32', mode='w+', order = 'F', shape=np.shape(beta_map))
	mm[:] = beta_map[:]


