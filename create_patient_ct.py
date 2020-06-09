import ROOT
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from multiprocessing import Pool
import glob
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator

number_format = "float"
scanner = 'ct_interpolation_test'
bytes_per_pixel = "4"

patient_ct_binary = '/home/baran/git/Simulations_GATE/CT_data/cirs_modified_for_ac.raw'
patient_ct_header = '/home/baran/git/Simulations_GATE/CT_data/cirs_modified.mhd'

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 
fov_z_dualhead_3x4 = 400. 

# CNAO FOV size in mm (needed to calculate the matrix size)
cnao_fov_x = 99.2 
cnao_fov_y = 220.8 
cnao_fov_z = 249.6 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 2.5
voxel_y = 2.5
voxel_z = 2.5 

# CNAO VOXEL size in mm (scalling factor: mm/pixel)
cnao_voxel_x = 1.6 
cnao_voxel_y = 1.6 
cnao_voxel_z = 1.6 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)
mat_z_dualhead_3x4 = int(fov_z_dualhead_3x4/voxel_z)

# CNAO MATRIX size 
cnao_mat_x = int(cnao_fov_x/cnao_voxel_x) 
cnao_mat_y = int(cnao_fov_y/cnao_voxel_y)
cnao_mat_z = int(cnao_fov_z/cnao_voxel_z)

#CT DATA
ct_voxel_x = 0.9765625 
ct_voxel_y = 0.9765625 
ct_voxel_z = 1.2

ct_mat_x = 512 
ct_mat_y = 512
ct_mat_z = 202

isocenter_x = -20.0 
isocenter_y = -230.0
isocenter_z = -390.0

offset_x = -249.512
offset_y = -478.512
offset_z = -530.2


# ########################
# PHANTOM TO FILE
# ########################

def load_and_convert_ct_to_lac():

	print 'Loading CT ...'
	ct = np.fromfile(patient_ct_binary, dtype='float32').reshape((ct_mat_x, ct_mat_y, ct_mat_z), order = 'F')
	print '... done!'
	print 'CT to umap transformation ...'
	ct_offset = -1000
#	ct = ct - ct_offset
	umap = np.zeros((ct.shape))
	#print umap.shape	
	#print np.min(ct)
	#print np.max(ct)
	# CT voltage tube 120kVp
	# All factors are taken from the Carney et al., Transforming CT images for attenuation correction in PET/CT, Medical Physics, Vol. 33, No. 4, 2006
	# Below the breaking point: u = below_a*(HU+1000)
	# Above the breaking point: u = above_a*(HU+1000) + above_b
	# Check if loaded CT is already from 0 HU or from -1000 HU
	break_point = 47	
#	break_point = 1047	
	below_a = 9.6E-5
	below_b = 0.0
	above_a = 5.10E-5
	above_b = 4.71E-2
	
	umap1 = below_a*(ct-ct_offset)+below_b
	umap2 = above_a*(ct-ct_offset)+above_b
#	umap1 = below_a*ct+below_b
#	umap2 = above_a*ct+above_b
	umap[ct<break_point] = umap1[ct<break_point]
	umap[ct>=break_point] = umap2[ct>=break_point]
	print '... done!'
	#print np.min(umap)
	#print np.max(umap)

	return umap


def save_umap_to_binary (umap):

	# CALCULATING ORIGINAL CT COORDINATES
	x=np.linspace( offset_x-isocenter_x , offset_x-isocenter_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
	y=np.linspace( offset_y-isocenter_y , offset_y-isocenter_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
	z=np.linspace( offset_z-isocenter_z , offset_z-isocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)
	'''
	print x
	print y 
	print z
	'''	
	linear_interp = RegularGridInterpolator((x, y, z), umap, bounds_error=False, fill_value=0.0)

# BARRELS, DUALHEAD_1x6, DUALHEAD_2x6 SETUPS

	print 'Processing umap for barrels and two dualheads ...'
	# CALCULATING NEW umap COORDINATES

	pts_normal = np.zeros((mat_x * mat_y * mat_z, 3))
	rescaled_umap_normal = np.zeros((mat_x, mat_y, mat_z))
	
	for i in range (0, rescaled_umap_normal.shape[2]):
		for j in range (0, rescaled_umap_normal.shape[1]):
			for k in range (0, rescaled_umap_normal.shape[0]):

					ind =   k   +   j*rescaled_umap_normal.shape[0]   +  i*rescaled_umap_normal.shape[0]*rescaled_umap_normal.shape[1]
					pts_normal[ind,0] = -fov_x/2 + voxel_x*(1 + k) - voxel_x/2
					pts_normal[ind,1] = -fov_y/2 + voxel_y*(1 + j) - voxel_y/2
					pts_normal[ind,2] = -fov_z/2 + voxel_z*(1 + i) - voxel_z/2

	#print pts_normal
	# TRILINEAR INTERPOLATION
	rescaled_umap = linear_interp(pts_normal).reshape(rescaled_umap_normal.shape, order = 'F')

	# SAVING NEW umap
	rescaled_umap.flatten('F').astype('float32').tofile('rescaled_cirs_umap.raw')
	print '... done!'

# DUALHEAD_3x4 SETUP

	print 'Processing umap for dualhead_3x4 ...'
	# CALCULATING NEW umap COORDINATES

	pts_normal = np.zeros((mat_x * mat_y * mat_z_dualhead_3x4, 3))
	rescaled_umap_dualhead_3x4 = np.zeros((mat_x, mat_y, mat_z_dualhead_3x4))
	
	for i in range (0, rescaled_umap_dualhead_3x4.shape[2]):
		for j in range (0, rescaled_umap_dualhead_3x4.shape[1]):
			for k in range (0, rescaled_umap_dualhead_3x4.shape[0]):

					ind =   k   +   j*rescaled_umap_dualhead_3x4.shape[0]   +  i*rescaled_umap_dualhead_3x4.shape[0]*rescaled_umap_dualhead_3x4.shape[1]
					pts_normal[ind,0] = -fov_x/2 + voxel_x*(1 + k) - voxel_x/2
					pts_normal[ind,1] = -fov_y/2 + voxel_y*(1 + j) - voxel_y/2
					pts_normal[ind,2] = -fov_z_dualhead_3x4/2 + voxel_z*(1 + i) - voxel_z/2

	#print pts_normal
	# TRILINEAR INTERPOLATION
	rescaled_umap = linear_interp(pts_normal).reshape(rescaled_umap_dualhead_3x4.shape, order = 'F')

	# SAVING NEW umap
	rescaled_umap.flatten('F').astype('float32').tofile('rescaled_cirs_umap_dualhead_3x4.raw')
	print '... done!'

# CNAO SETUP

	print 'Processing umap for cnao_lso ...'
	pts_cnao = np.zeros((cnao_mat_x * cnao_mat_y * cnao_mat_z, 3))
	cnao_umap = np.zeros((cnao_mat_x, cnao_mat_y, cnao_mat_z))
	#print cnao_pmma.shape

	
	for i in range (0, cnao_umap.shape[2]):
		for j in range (0, cnao_umap.shape[1]):
			for k in range (0, cnao_umap.shape[0]):

					ind =   k   +   j*cnao_umap.shape[0]   +  i*cnao_umap.shape[0]*cnao_umap.shape[1]
					pts_cnao[ind,0] = -cnao_fov_x/2 + (cnao_voxel_x)*(1 + k) - (cnao_voxel_x/2)
					pts_cnao[ind,1] = -cnao_fov_y/2 + (cnao_voxel_y)*(1 + j) - (cnao_voxel_y/2)
					pts_cnao[ind,2] = -cnao_fov_z/2 + (cnao_voxel_z)*(1 + i) - (cnao_voxel_z/2) 

	
	# TRILINEAR INTERPOLATION
	rescaled_umap = linear_interp(pts_cnao).reshape(cnao_umap.shape, order = 'F')

	# SAVING NEW umap
	rescaled_umap.flatten('F').astype('float32').tofile('rescaled_cirs_umap_cnao_lso.raw')

	print '... done!'




def save_rescaled_umap_to_header ():

        with open('rescaled_umap.hdr', "w+") as f:
		f.write("!name of data  file :=  {0}\n".format("rescaled_cirs_umap.raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("imagedata  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(mat_x))
		f.write("!matrix  size  [2] := {0}\n".format(mat_y))
		f.write("!matrix  size  [3] := {0}\n".format(mat_z))
		f.write("!number  format  := float\n")
		f.write("!number  of bytes  per  pixel  := {0}\n".format(bytes_per_pixel))
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(voxel_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(voxel_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(voxel_z))
		f.write("image  duration (sec) := 1\n")

def save_rescaled_umap_dualhead_3x4_to_header ():

        with open('rescaled_umap_dualhead_3x4.hdr', "w+") as f:
		f.write("!name of data  file :=  {0}\n".format("rescaled_cirs_umap_dualhead_3x4.raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("imagedata  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(mat_x))
		f.write("!matrix  size  [2] := {0}\n".format(mat_y))
		f.write("!matrix  size  [3] := {0}\n".format(mat_z_dualhead_3x4))
		f.write("!number  format  := float\n")
		f.write("!number  of bytes  per  pixel  := {0}\n".format(bytes_per_pixel))
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(voxel_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(voxel_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(voxel_z))
		f.write("image  duration (sec) := 1\n")
	
def save_cnao_umap_to_header ():

        with open('rescaled_umap_cnao_lso.hdr', "w+") as f:
		f.write("!name of data  file :=  {0}\n".format("rescaled_cirs_umap_cnao_lso.raw"))
		f.write("!total  number  of  images  := 1\n")
		f.write("imagedata  byte  order :=  LITTLEENDIAN\n")
		f.write("number  of  dimensions  := 3\n")
		f.write("!matrix  size  [1] := {0}\n".format(cnao_mat_x))
		f.write("!matrix  size  [2] := {0}\n".format(cnao_mat_y))
		f.write("!matrix  size  [3] := {0}\n".format(cnao_mat_z))
		f.write("!number  format  := float\n")
		f.write("!number  of bytes  per  pixel  := {0}\n".format(bytes_per_pixel))
		f.write("scaling  factor (mm/pixel) [1] := {0}\n".format(cnao_voxel_x))
		f.write("scaling  factor (mm/pixel) [2] := {0}\n".format(cnao_voxel_y))
		f.write("scaling  factor (mm/pixel) [3] := {0}\n".format(cnao_voxel_z))
		f.write("image  duration (sec) := 1\n")



# ########################
# MAIN
# ########################


if __name__ == '__main__':


	umap = load_and_convert_ct_to_lac()
	save_umap_to_binary(umap)
	save_rescaled_umap_to_header()
	save_rescaled_umap_dualhead_3x4_to_header()
#	save_cnao_umap_to_header()
	

