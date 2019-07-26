'''
Created on 5 Sep. 2017

@author: jakubb
'''

import ROOT
import sys
import os
#import numpy as np
#import matplotlib.pyplot as plt
#from multiprocessing import Pool
#import glob

#back_to_back_number = 100000000000
#number_format = "float"
#bytes_per_pixel = "4"
#scanner = "barrel"
#phantom_name = "pmma_5_5_20_cm_shifted"

#sensitivity_binary_castor = "/home/baran/git/CastorTools/castor_sensitivity_maps/sensitivity_"+scanner+".raw"
#sensitivity_header_castor = "/home/baran/git/CastorTools/castor_sensitivity_maps/sensitivity_"+scanner+".hdr"

# FOV size in mm (needed to calculate the matrix size)
#fov_x = 400. 
#fov_y = 400. 
#fov_z = 500. 

# VOXEL size in mm (scalling factor: mm/pixel)
#voxel_x = 5. 
#voxel_y = 5. 
#voxel_z = 5. 

# MATRIX size 
#mat_x = int(fov_x/voxel_x) 
#mat_y = int(fov_y/voxel_y)
#mat_z = int(fov_z/voxel_z)

# Annihilation number per voxel
#prim_per_voxel = float(back_to_back_number)/(mat_x*mat_y*mat_z)

# Maps path
#sens_maps_path= '/home/baran/git/CastorTools/castor_sensitivity_maps/barrel_data'


# ########################
# Y/Z COORDINATE SMOOTHING
# ########################
def smooth (val, randomizer, smoothing_sigma):

	new_val = randomizer (val, smoothing_sigma)
	if new_val <= -250.:
		new_val = -249.9
	if new_val >= 250.:
		new_val = 249.9

	return new_val




# ########################
# MAIN
# ########################


if __name__ == '__main__':

	inbeam = 1
	smoothing = 0
	barrel = 1	# 1 correspond to the BARREL setup; 0 - DUALHEAD
	smoothing_sigma = 11 # in mm  

	discretization_dim = 5. # in mm; correspond to the voxel dimension in the 
	plastic_len = 500. # in mm;
	
	time_gap = 60. # in seconds: sum of the irradiation time and time needed to transport the patient 

	root_dir_path = "/home/baran/Desktop/castor_recons/barrel_cubic_phantom/"
	file_name = "cubic_phantom_results.root"
	root_file_path = root_dir_path + file_name
	new_root_file_path = root_dir_path + "inbeam.root"
	if not inbeam:
		new_root_file_path = root_dir_path + "offbeam.root"

	try:
		old = ROOT.TFile.Open(root_file_path)
		new = ROOT.TFile(new_root_file_path, "RECREATE")
	except:
		print "Error opening the files."
		sys.exit(2)

	old_tree = old.Get("Coincidences")
	new_tree = old_tree.CloneTree(0)

	randomizer = ROOT.gRandom.Gaus

	for i, entry in enumerate(old_tree):

		if (smoothing == 1):
			if (barrel):
				entry.globalPosZ1 = smooth(entry.globalPosZ1, randomizer, smoothing_sigma)
				entry.globalPosZ2 = smooth(entry.globalPosZ2, randomizer, smoothing_sigma)
				entry.layerID1 = int ( (entry.globalPosZ1+(plastic_len/2.)) / discretization_dim)
				entry.layerID2 = int ( (entry.globalPosZ2+(plastic_len/2.)) / discretization_dim)
			else:
				entry.globalPosY1 = smooth(entry.globalPosY1, randomizer, smoothing_sigma)
				entry.globalPosY2 = smooth(entry.globalPosY2, randomizer, smoothing_sigma)
				entry.layerID1 = int ( (entry.globalPosY1+(plastic_len/2.)) / discretization_dim)
				entry.layerID2 = int ( (entry.globalPosY2+(plastic_len/2.)) / discretization_dim)

			if (entry.layerID1 == 100):
				entry.layerID1 = 99
			if (entry.layerID2 == 100):
				entry.layerID2 = 99

		if (inbeam == 1 and entry.time1 < time_gap and entry.time2 < time_gap):
			new_tree.Fill()
		if (inbeam == 0 and entry.time1 >= time_gap and entry.time2 >= time_gap):
			new_tree.Fill()
		#if (entry.crystalID1 == 10):
		#	print "CrystalID1: {0}, layerID1: {1}, globalPosZ1: {2}\n".format(entry.crystalID1, entry.layerID1, entry.globalPosZ1)
			 
			
	
	
	#new_tree.Print()
	new_tree.Write()
	new.Save()
	old.Close()
	new.Close()





