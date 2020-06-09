import ROOT
import sys
import os
from multiprocessing import Pool
from subprocess import check_call as p
import numpy as np
import glob


sensitivity_dir_path = '/home/baran/git/CastorTools/castor_sensitivity_maps_data/'
scanner_name = 'barrel2'
dir_path = sensitivity_dir_path + scanner_name + "_data/" 
cpu_pool =  20

def run_shell (sensitivity_file):

	'''
	Run the command from the prm code dedicated to produce additional LORs (code available here: https://github.com/jawka/GateTools)
	'''

	shell_command = 'prm --dir_path ' + dir_path + ' --scanner_name ' + scanner_name + ' --mode 12 --res_blur_mode '  + scanner_name + ' --coincidences_file ' + sensitivity_file 
	#print shell_command
	p(shell_command, bufsize=0, shell=True)


# ########################
# MAIN
# ########################


if __name__ == '__main__':


	'''

	Script to copy and produce additional LORs to cover whole FOV 
	in order to produce merged sensitivity and attenuation map	

	'''

	os.chdir(sensitivity_dir_path+scanner_name+'_data/')
	lm_list = glob.glob('*.root')
	print len(lm_list)

	p = Pool(cpu_pool)
	p.map (run_shell, lm_list)
	#run_shell (lm_list)
	print "Done!"

