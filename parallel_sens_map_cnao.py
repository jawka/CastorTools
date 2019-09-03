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


back_to_back_number = 100000000000
number_format = 'float'
bytes_per_pixel = '4'
scanner = 'cnao_lso'

sensitivity_binary_castor = "/home/baran/git/CastorTools/castor_sensitivity_maps_data/sensitivity_"+scanner+".raw"
sensitivity_header_castor = "/home/baran/git/CastorTools/castor_sensitivity_maps_data/sensitivity_"+scanner+".hdr"

# FOV size in mm (needed to calculate the matrix size)
#fov_x = 400. 
#fov_y = 400. 
#fov_z = 400. 
# dualhead 3x4
#fov_x = 400. 
#fov_y = 400. 
#fov_z = 400. 
# cnao_lso
fov_x = 99.2 
fov_y = 220.8 
fov_z = 249.6 


# VOXEL size in mm (scalling factor: mm/pixel)
#voxel_x = 5. 
#voxel_y = 5. 
#voxel_z = 5. 

# cnao_lso
voxel_x = 1.6 
voxel_y = 1.6 
voxel_z = 1.6 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)

# Annihilation number per voxel
prim_per_voxel = float(back_to_back_number)/(mat_x*mat_y*mat_z)

# Maps path
sens_maps_path= '/home/baran/git/CastorTools/castor_sensitivity_maps_data/'+scanner+'_data'


# ########################
# SENSITIVITY TO FILE
# ########################

def save_sensitivity_to_binary_argument (sensitivity):

	print 'Bytes order: {0} (should be \'little\'). If \'big\', bytes swapping will be automatically applied.'.format(sys.byteorder)
	if sys.byteorder != 'little':
		print 'Byte swapping ...'
		sensitivity.byteswap()
		print 'DONE'

	#np.save(sensitivity_binary_castor, sensitivity)
	sensitivity.flatten('F').astype('float32').tofile(sensitivity_binary_castor)

def save_sensitivity_to_header ():

        with open(sensitivity_header_castor, "w+") as f:
		f.write("!name of data  file :=  {0}\n".format("sensitivity_"+scanner+".raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("imagedata  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(sensitivity.shape[0]))
		f.write("!matrix  size  [2] := {0}\n".format(sensitivity.shape[1]))
		f.write("!matrix  size  [3] := {0}\n".format(sensitivity.shape[2]))
		f.write("!number  format  := "+number_format+"\n")
		f.write("!number  of bytes  per  pixel  := {0}\n".format(bytes_per_pixel))
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(voxel_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(voxel_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(voxel_z))
		f.write("image  duration (sec) := 1\n")

def create_sensitivity_parallel (results_path):


	batch_number = results_path.split('/')[-1].split('results')[-1].split('.')[0]

	sensitivity = np.zeros((mat_x, mat_y, mat_z))

	try:
		f = ROOT.TFile.Open(results_path)
		print "Processed map: {0}".format(batch_number)
		#print "File {0} opened successfully!\n".format(results_path)
	except:
		print "File {0} doesn't open properly or doesn't exist!".format(results_path)
		sys.exit(2)
 	
	tree = f.Get("Coincidences")
	entries = tree.GetEntries()

#	print "Number of coincidences: {0}".format(entries) 

	temp_x = 0.	
	temp_y = 0.
	temp_z = 0.
	for i, entry in enumerate(tree):
	#	if ((i+1)%100000 == 0):
	#		print "{0}00 k events analysed".format(int((i+1)/100000))
		if (entry.sourcePosX1 == entry.sourcePosX2 and entry.sourcePosY1 == entry.sourcePosY2 and entry.sourcePosZ1 == entry.sourcePosZ2 and entry.eventID1 == entry.eventID2):
			#print "Annihilation position (in mm): {0} {1} {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosX1)
			temp_x = int( (entry.sourcePosX1 + fov_x/2) / voxel_x)
			temp_y = int( (entry.sourcePosY1 + fov_y/2) / voxel_y)
			temp_z = int( (entry.sourcePosZ1 + fov_z/2) / voxel_z)
			if (entry.sourcePosX1 == fov_x/2.):
				print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
				temp_x -= 1
			if (entry.sourcePosY1 == fov_y/2.):
				print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
				temp_y -= 1
			if (entry.sourcePosZ1 == fov_z/2.):
				print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
				temp_z -= 1
			#print "Corresponding matrix coordinates: {0} {1} {2}\n".format( temp_x, temp_y, temp_z)
			#if (temp_x >= (mat_x/2) or temp_y >= (mat_y/2) or temp_z >= (mat_z/2)):
			#	print "Source position: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1) 
			sensitivity[temp_x, temp_y, temp_z] += 1

	f.Close()
	plot_flag = 0
	if plot_flag == 1:
		fig, axs = plt.subplots(2,2)
	
		im1 = axs[0,0].imshow(sensitivity[31,:,:])
		axs[0,0].set_title("x == 31")
		fig.colorbar(im1, axs[0,0])

		im2 = axs[0,1].imshow(sensitivity[:,69,:])
		axs[0,1].set_title("y == 69")
		fig.colorbar(im2, axs[0,1])
	
		im3 = axs[1,0].imshow(sensitivity[:,:,78])
		axs[1,0].set_title("z == 78")
		fig.colorbar(im3, axs[1,0])
	
#		im4 = axs[1,1].imshow(sensitivity[:,:,75])
#		axs[1,1].set_title("75")
#		fig.colorbar(im4, axs[1,1])
	
		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		fig.colorbar(im1, cax=cbar_ax)
	
		plt.show()


	sensitivity_file = sens_maps_path+'/sensitivity_maps/sens_map_{0}'.format(batch_number)	
	np.save(sensitivity_file, sensitivity)


# ########################
# MAIN
# ########################


if __name__ == '__main__':




	print "FOV size: {0} x {1} x {2} mm3".format(fov_x, fov_y, fov_z)
	print "Voxel size: {0} x {1} x {2} (mm/pixel)3".format(voxel_x, voxel_y, voxel_z)
	print "Matrix size: {0} x {1} x {2}\n".format(mat_x, mat_y, mat_z)

	# Sensitivity matrix
	sensitivity = np.zeros((mat_x, mat_y, mat_z))

	print "Full sensitivity matrix size: {0}".format(sensitivity.shape)
	print "Simulated coincidences number: {0} G".format(back_to_back_number/1000000000)
	print "Simulated coincidences per voxel number: {0}\n".format(prim_per_voxel)
	

	## SENSITIVITY MAPS PROCESSING


	maps_list = glob.glob(sens_maps_path+'/*.root')
#	maps_list = glob.glob(sens_maps_path+'/*.root')[0]   # number 1917
	print len(maps_list)
	print maps_list	
	if not os.path.exists(sens_maps_path+'/sensitivity_maps'):
		os.mkdir (sens_maps_path+'/sensitivity_maps', 0775)
	p = Pool(35)
	p.map (create_sensitivity_parallel, maps_list)
	#create_sensitivity_parallel (maps_list)
	matrix_list = glob.glob(sens_maps_path+'/sensitivity_maps/sens_map_*')


	print "FOV size: {0} x {1} x {2} mm3".format(fov_x, fov_y, fov_z)
	print "Voxel size: {0} x {1} x {2} (mm/pixel)3".format(voxel_x, voxel_y, voxel_z)
	print "Matrix size: {0} x {1} x {2}\n".format(mat_x, mat_y, mat_z)


	print 'Number of matrices to be analysed: {}'.format(len(matrix_list))
	
	for i, val in enumerate (matrix_list):

		batch_number = matrix_list[i].split('/')[-1].split('results')[-1].split('.')[0]
#		print "Loaded matrix: {0}".format(batch_number)
		temp_matrix = np.load(matrix_list[i])
		sensitivity = np.add(sensitivity, temp_matrix)

	print "\nSensitivity matrix counts:"	
	print "Sensitivity matrix min: {}".format(sensitivity.min())
	print "Sensitivity matrix max: {}\n".format(sensitivity.max())

	sensitivity[sensitivity < 1] = 1

	sensitivity = sensitivity/prim_per_voxel
	os.chdir(sens_maps_path)
	sensitivity_path = 'final_sensitivity_map'

	print "Sensitivity matrix min: {}".format(sensitivity.min())
	print "Sensitivity matrix max: {}\n".format(sensitivity.max())

	print "\nSensitivity matrix counts/voxel:"	
	print "Sensitivity matrix min: {}".format(sensitivity.min())
	print "Sensitivity matrix max: {}\n".format(sensitivity.max())
	np.save(sensitivity_path, sensitivity)
	

#	fig, axs = plt.subplots(2,2)

#	im1 = axs[0,0].imshow(sensitivity[:,:,24])
#	axs[0,0].set_title("24")
#	fig.colorbar(im1, axs[0,0])

#	im2 = axs[0,1].imshow(sensitivity[:,:,73])
#	axs[0,1].set_title("73")
#	fig.colorbar(im2, axs[0,1])

#	im3 = axs[1,0].imshow(sensitivity[:,:,74])
#	axs[1,0].set_title("74")
#	fig.colorbar(im3, axs[1,0])

#	im4 = axs[1,1].imshow(sensitivity[:,:,75])
#	axs[1,1].set_title("75")
#	fig.colorbar(im4, axs[1,1])

#	fig.subplots_adjust(right=0.8)
#	cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#	fig.colorbar(im1, cax=cbar_ax)

#	plt.show()



	save_sensitivity_to_binary_argument(sensitivity)
	save_sensitivity_to_header()



