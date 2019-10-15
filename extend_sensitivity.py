import ROOT
import sys
import os
from multiprocessing import Pool
from subprocess import check_call as p
import numpy as np
import glob


sensitivity_dir_path = '/home/baran/git/CastorTools/castor_sensitivity_maps_data/'
scanner_name = 'dualhead_2x6'
dir_path = sensitivity_dir_path + scanner_name + "_data/" 
cpu_pool =  35
extended_sensitivity_list_file = '/home/baran/Desktop/castor_recons/barrel_first_voxel_5_5_5/extended_sensitivity_list.txt'

def run_shell (sensitivity_file):

	shell_command = 'prm --dir_path ' + dir_path + ' --scanner_name ' + scanner_name + ' --mode 12 --res_blur_mode '  + scanner_name + ' --coincidences_file ' + sensitivity_file 
	#print shell_command
	p(shell_command, bufsize=0, shell=True)


# ########################
# MAIN
# ########################


if __name__ == '__main__':


	os.chdir(sensitivity_dir_path+scanner_name+'_data/')
	lm_list = glob.glob('*.root')
	print len(lm_list)

	p = Pool(cpu_pool)
	p.map (run_shell, lm_list)
	#run_shell (lm_list)
	print "Done!"


#	lm_list = glob.glob(sensitivity_dir_path+scanner_name+'_data/extended*.root')
#	print len(lm_list)	
#	with open (extended_sensitivity_list_file, 'w') as f:
#		for i, ll in enumerate (lm_list):
#			f.write(ll+'\n')

#	prm --dir_path ~/git/CastorTools/castor_sensitivity_maps_data/test_extend_dualhead_1x6/ --scanner_name dualhead_1x6 --mode 12 --res_blur_mode dualhead_1x6 --coincidences_file sensitivity_results100.root
