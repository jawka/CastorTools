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

number_format = "float"
bytes_per_pixel = "4"
scanner = "dualhead"
phantom_name = "test_1_8"


phantom_binary_castor = "/home/baran/git/CastorTools/castor_phantoms/"+phantom_name+".raw"
phantom_header_castor = "/home/baran/git/CastorTools/castor_phantoms/"+phantom_name+".hdr"



# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 5. 
voxel_y = 5. 
voxel_z = 5. 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)


# ########################
# PHANTOM TO FILE
# ########################

def save_phantom_to_binary ():

	phantom = np.zeros((mat_x, mat_y, mat_z))

	pmma_HU = 0.10473 # linear attenation coefficient for PMMA in cm-1
	pmma_phantom = np.full((10,10,40), pmma_HU)
	#print pmma_phantom.shape
	#phantom [35:45, 35:45, 30:70] = pmma_phantom
		

	water_HU = 0.1 # linear attenation coefficient for PMMA in cm-1
	pmma_phantom = np.full((50,50,60), water_HU)
	#print pmma_phantom.shape
	#phantom [15:65, 15:65, 20:80] = pmma_phantom

	test_value = 0.00001 # linear attenation coefficient for PMMA in cm-1
	test_1_8_phantom = np.full((40,40,50), test_value)
	#print pmma_phantom.shape
	phantom [40:, 40:, 50:] = test_1_8_phantom

	print 'Bytes order: {0} (should be \'little\'). If \'big\', bytes swapping will be automatically applied.'.format(sys.byteorder)
	if sys.byteorder != 'little':
		print 'Byte swapping ...'
		sensitivity.byteswap()
		print 'DONE'

	#np.save(sensitivity_binary_castor, sensitivity)
	#pmma_phantom.flatten('F').astype('float32').tofile(phantom_binary_castor)

	#phantom = np.ones((mat_x, mat_y, mat_z))
	phantom.flatten('F').astype('float32').tofile(phantom_binary_castor)


def save_phantom_to_header ():

        with open(phantom_header_castor, "w+") as f:
		f.write("!name of data  file :=  {0}\n".format(phantom_name+".raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("imagedata  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(mat_x))
		f.write("!matrix  size  [2] := {0}\n".format(mat_y))
		f.write("!matrix  size  [3] := {0}\n".format(mat_z))
		f.write("!number  format  := "+number_format+"\n")
		f.write("!number  of bytes  per  pixel  := {0}\n".format(bytes_per_pixel))
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(voxel_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(voxel_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(voxel_z))
		f.write("image  duration (sec) := 1\n")





# ########################
# MAIN
# ########################


if __name__ == '__main__':




	save_phantom_to_binary()
	save_phantom_to_header()


