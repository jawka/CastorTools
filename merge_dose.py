import sys
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

def read_dose_map_hdr(hdr_file):

	f1 = open('dose-Dose.mhd', 'w')
	matrix_size = np.zeros((3,1))
	voxel_size = np.zeros((3,1))
	with open(hdr_file, 'r')as f:
		for line in f:
			if line.startswith('ElementSpacing'):
				voxel_size = map(float, line.split(' = ')[-1].split(" "))
			if line.startswith('DimSize'):
				matrix_size = map(int, line.split(' = ')[-1].split(" "))
			if line.startswith('ElementDataFile'):
				f1.write('ElementDataFile = dose-Dose.raw\n')
			else:
				f1.write(line)

	f1.close()
	print matrix_size
	print voxel_size
	return (matrix_size, voxel_size)


if __name__ == '__main__':
	
	'''

	Script to merge dose maps from splitted simulation

	'''

	dose_dir_path = '/home/baran/git/Simulations_GATE/patient_positrons'
	
	os.chdir(dose_dir_path)
	dose_hdr = glob.glob('*Dose.mhd')
	dose_raw = glob.glob('*Dose.raw')
	print dose_hdr
	matrix_size, voxel_size = read_dose_map_hdr(dose_hdr[0])
	sum_dose = np.zeros((matrix_size))

	for ff in range(1,21):
		if os.path.exists('dose'+str(ff)+'-Dose.raw'):
			print ff
			temp_dose = np.fromfile('dose'+str(ff)+'-Dose.raw', dtype = 'float32').reshape(matrix_size)
			sum_dose = sum_dose+temp_dose
		else:
			print '{} dose file is missing!'.format (ff)
	
	sum_dose.flatten().astype('float32').tofile('dose-Dose.raw')
	
	

