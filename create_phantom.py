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
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator


number_format = "float"
bytes_per_pixel = "4"
scanner = "dualhead"
phantom_name = "cnao_water_phantom"


phantom_binary_castor = "/home/baran/git/CastorTools/castor_phantoms/"+phantom_name+".raw"
phantom_header_castor = "/home/baran/git/CastorTools/castor_phantoms/"+phantom_name+".hdr"



# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 

# CNAO FOV size in mm (needed to calculate the matrix size)
cnao_fov_x = 99.2 
cnao_fov_y = 220.8 
cnao_fov_z = 249.6 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 5. 
voxel_y = 5. 
voxel_z = 5. 

# CNAO VOXEL size in mm (scalling factor: mm/pixel)
cnao_voxel_x = 1.6 
cnao_voxel_y = 1.6 
cnao_voxel_z = 1.6 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)


# CNAO MATRIX size 
cnao_mat_x = int(cnao_fov_x/cnao_voxel_x) 
cnao_mat_y = int(cnao_fov_y/cnao_voxel_y)
cnao_mat_z = int(cnao_fov_z/cnao_voxel_z)


# ########################
# PHANTOM TO FILE
# ########################

def save_phantom_to_binary ():

	phantom = np.zeros((mat_x, mat_y, mat_z))

	## PMMA PHANTOM

	pmma_HU = 0.10473 # linear attenation coefficient for PMMA in cm-1
	pmma_phantom = np.full((10,10,40), pmma_HU)
	#print pmma_phantom.shape

	phantom [35:45, 35:45, 30:70] = pmma_phantom

	#FOR DUALHEAD_3x4 ONLY
	#phantom [35:45, 35:45, 20:60] = pmma_phantom







	#FOR CNAO INTERPOLATIOM

	#matrix_5mm = np.load("5mm_cubic_voxel.npy")
	#matrix_10mm = np.load("10mm_cubic_voxel.npy")

	#print matrix_10mm.shape

	x=np.linspace(-fov_x/2+(voxel_x/2),fov_x/2-(voxel_x/2),phantom.shape[0])
	y=np.linspace(-fov_y/2+(voxel_y/2),fov_y/2-(voxel_y/2),phantom.shape[1])
	z=np.linspace(-fov_z/2+(voxel_z/2),fov_z/2-(voxel_z/2),phantom.shape[2])

	#x[0] = -fov_x/2
	#y[0] = -fov_y/2
	#z[0] = -250.

	#x[-1] = fov_x/2
	#y[-1] = fov_y/2
	#z[-1] = 250.

	#print x
	#print y 
	#print z
	
	linear_interp = RegularGridInterpolator((x, y, z), phantom)
	#print phantom.shape

	pts = np.zeros((cnao_mat_x * cnao_mat_y * cnao_mat_z, 3))
	#print pts.shape

	cnao_pmma = np.zeros((cnao_mat_x, cnao_mat_y, cnao_mat_z))
	#print cnao_pmma.shape

	
	for i in range (0, cnao_pmma.shape[2]):
		for j in range (0, cnao_pmma.shape[1]):
			for k in range (0, cnao_pmma.shape[0]):

					ind =   k   +   j*cnao_pmma.shape[0]   +  i*cnao_pmma.shape[0]*cnao_pmma.shape[1]
# FOR CNAO interpolation
					pts[ind,0] = -cnao_fov_x/2 + (cnao_voxel_x)*(1 + k) - (cnao_voxel_x/2)
					pts[ind,1] = -cnao_fov_y/2 + (cnao_voxel_y)*(1 + j) - (cnao_voxel_y/2)
					pts[ind,2] = -cnao_fov_z/2 + (cnao_voxel_z)*(1 + i) - (cnao_voxel_z/2) 

	#print pts

	cnao_pmma = linear_interp(pts).reshape(cnao_pmma.shape, order = 'F')
	#phantom = cnao_pmma
	













		

	## WATER PHANTOM

	water_HU = 0.1 # linear attenation coefficient for PMMA in cm-1
	#water_phantom = np.full((50,50,60), water_HU)
	#print water_phantom.shape
	
	#phantom [15:65, 15:65, 20:80] = water_phantom
	#FOR DUALHEAD_3x4 ONLY
	#phantom [15:65, 15:65, 10:70] = water_phantom
	#FOR CNAO
	phantom = np.zeros((cnao_mat_x, cnao_mat_y, cnao_mat_z))
	water_phantom = np.full((40,100,100), water_HU)
	phantom [11:51, 19:119, 28:128] = water_phantom


	## TESTING

	test_value = 0.00001 # linear attenation coefficient for PMMA in cm-1
	test_1_8_phantom = np.full((40,40,50), test_value)
	#print pmma_phantom.shape
	#phantom [40:, 40:, 50:] = test_1_8_phantom


	## SAVING

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

	
def save_cnao_phantom_to_header ():

        with open(phantom_header_castor, "w+") as f:
		f.write("!name of data  file :=  {0}\n".format(phantom_name+".raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("imagedata  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(cnao_mat_x))
		f.write("!matrix  size  [2] := {0}\n".format(cnao_mat_y))
		f.write("!matrix  size  [3] := {0}\n".format(cnao_mat_z))
		f.write("!number  format  := "+number_format+"\n")
		f.write("!number  of bytes  per  pixel  := {0}\n".format(bytes_per_pixel))
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(cnao_voxel_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(cnao_voxel_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(cnao_voxel_z))
		f.write("image  duration (sec) := 1\n")



# ########################
# MAIN
# ########################


if __name__ == '__main__':




	save_phantom_to_binary()
	#save_phantom_to_header()
	save_cnao_phantom_to_header()
	

