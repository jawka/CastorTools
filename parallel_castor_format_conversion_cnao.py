import ROOT
import sys
import os
from multiprocessing import Pool
from subprocess import check_call as p
import numpy as np
import glob


castor_scanner = 'CNAO_LSO'
gate_scanner = 'cnao_lso'
sensitivity_dir_path = '/home/baran/git/CastorTools/castor_sensitivity_maps_data/' + gate_scanner + '_data/'

tof_flag = 0

cpu_pool =  5

def run_shell (sensitivity_file):

	no = sensitivity_file.split('results')[-1].split('.')[0]
	shell_command = 'castor-GATERootToCastor -t -i ' + sensitivity_file + ' -s ' + castor_scanner + ' -vb 0 -m ~/git/Simulations_GATE/system_matrix_' + gate_scanner + '/' + gate_scanner + '.mac' + ' -o extended_castor_' + no
	if tof_flag:
		shell_command = 'castor-GATERootToCastor -t -TOF_reso 500 -i ' + sensitivity_file + ' -s ' + castor_scanner + ' -vb 0 -m ~/git/Simulations_GATE/system_matrix_' + gate_scanner + '/' + gate_scanner + '.mac' + ' -o extended_castor_' + no
#	print shell_command
	p(shell_command, bufsize=0, shell=True)
	


# ########################
# MAIN
# ########################


if __name__ == '__main__':



	os.chdir(sensitivity_dir_path)
	lm_list = glob.glob('sensitivity_results*.root')

	## DATA CONVERSION

	print '\nNumber of list-mode files: {}'.format (len(lm_list))
	print 'Data conversion...' 
	pl = Pool(cpu_pool)
	pl.map (run_shell, lm_list)
	print "... done!\n"

	## HEADER CREATION
	
	all_lm_events = 0
	counter = 0
	duration = 0
	start_time = 0
	tof_reso = 0
	list_tof = 0
	
	print '\nList-mode events counting...' 
	if not os.path.exists('extended_castor_df.Cdh'):
		
		headers_list = glob.glob('extended_*.Cdh')
		print 'Number of header files: {}'.format (len(headers_list))
		for f_name in headers_list:
			with open(f_name, 'r') as hfile:
				lines=hfile.readlines()
#				line2=lines[1]
#				line5=lines[4]
#				line6=lines[5]
				all_lm_events = all_lm_events + int (lines[1].split(' ')[-1])
				start_time = int (lines[4].split(' ')[-1])
				duration = int (lines[5].split(' ')[-1])
				if tof_flag:
#					line11=lines[10]
#					line12=lines[11]
					tof_reso = float (lines[10].split(' ')[-1])
					list_tof_temp = float (lines[11].split(' ')[-1])
					if (list_tof_temp > list_tof):
						list_tof = list_tof_temp
				
		print 'All list-mode events: {0}'.format(all_lm_events)
		print "... done!"

		print 'Header creation...' 
		with open('extended_castor_df.Cdh', 'w') as hfile:

			hfile.write('Data filename: ' + gate_scanner + '_df.Cdf\n')
			hfile.write('Number of events: {0}\n'.format(all_lm_events))
			hfile.write('Data mode: list-mode\n')
			hfile.write('Data type: PET\n')
			hfile.write('Start time (s): {0}\n'.format(start_time))
			hfile.write('Duration (s): {0}\n'.format(duration))
			hfile.write('Scanner name: ' + castor_scanner+'\n')
			hfile.write('Calibration factor: 1\n')
			hfile.write('Isotope: unknown\n')
			if tof_flag:
				hfile.write('TOF information flag: 1\n')
				hfile.write('TOF resolution (ps): {0}\n'.format(tof_reso))
				hfile.write('List TOF measurement range (ps): {0}\n'.format(list_tof))


		print "... done!\n"
	else:
		print ' ... extended_castor_df.Cdh header file already exists !!!\n'

	## DATA CONCATENATION

	if not os.path.exists('extended_castor_df.Cdf'):
		print '\nData concatenation...' 

		shell_command = 'cat extended_castor_*_df.Cdf > extended_castor_df.Cdf'
		print shell_command
		p(shell_command, bufsize=0, shell=True)

		print "... done!\n"
	else:
		print ' ... extended_castor_df.Cdf list-mode file already exists !!!\n'

	


