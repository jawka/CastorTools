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
import scipy.ndimage.filters as img_filters


back_to_back_number = 1000000000
number_format = "float32"
bytes_per_pixel = "4"

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 

# FOV size in mm (needed to calculate the matrix size)
phantom_x = 50. 
phantom_y = 50. 
phantom_z = 200. 

# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 5. 
voxel_y = 5. 
voxel_z = 5. 

# VOXEL size in mm (scalling factor: mm/pixel)
dose_voxel_x = .5 
dose_voxel_y = .5 
dose_voxel_z = .5

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)

activity_dimx = int(phantom_x/voxel_x)
activity_dimy = int(phantom_y/voxel_y)
activity_dimz = int(phantom_z/voxel_z)


dose_dimx = int(phantom_x/dose_voxel_x)
dose_dimy = int(phantom_y/dose_voxel_y)
dose_dimz = int(phantom_z/dose_voxel_z)


dose_path_dir = '/home/baran/git/Simulations_GATE/proton_barrel/'
dose_name = 'dose_phantom' 

image_path = '/home/baran/Desktop/castor_recons/barrel_first_voxel_5_5_5/PET_images/PET_images_it2.img'

proton_beam_path='/home/baran/git/Simulations_GATE/proton_barrel/ps_actor_phantom'
positron_map_low_path='/home/baran/git/Simulations_GATE/proton_barrel/positron_map_low'
positron_map_high_path='/home/baran/git/Simulations_GATE/proton_barrel/positron_map_high'

# ########################
# CALCULATE THE DOSE
# ########################

def calculate_dose ():

	final_dose = np.zeros((dose_dimx, dose_dimy, dose_dimz))	
	for i in range (0, 20):
		image_path = dose_path_dir + dose_name + str(i+1) + '-Dose.raw'
		temp_image = np.fromfile(image_path, dtype=number_format).reshape((dose_dimx,dose_dimy,dose_dimz), order = 'F')
		#print temp_image.shape
		final_dose = final_dose + temp_image

	z_dose_high = final_dose.sum(axis=0).sum(axis=0)
	z_linspace_high = np.linspace(-phantom_z/2 + dose_voxel_z/2, phantom_z/2-dose_voxel_z/2, dose_dimz)


	z_dose_low = np.zeros((int(phantom_z/voxel_z)))
	for i in range(0, z_dose_high.shape[0]):
		ind = int(i/(voxel_z/dose_voxel_z))
		z_dose_low[ind] += z_dose_high[i]
	z_linspace_low = np.linspace(-phantom_z/2 + voxel_z/2, phantom_z/2-voxel_z/2, activity_dimz)

	plot_flag = 0
	if plot_flag == 1:

		fig, axs = plt.subplots(2,3)

		im1 = axs[0,0].imshow(final_dose[:,:,24])
		axs[0,0].set_title("24")
		fig.colorbar(im1, axs[0,0])
			
		im2 = axs[0,1].imshow(final_dose[:,:,73])
		axs[0,1].set_title("73")
		fig.colorbar(im2, axs[0,1])
	
		im3 = axs[1,0].imshow(final_dose[:,:,74])
		axs[1,0].set_title("74")
		fig.colorbar(im3, axs[1,0])
			
		im4 = axs[1,1].imshow(final_dose[:,:,75])
		axs[1,1].set_title("75")
		fig.colorbar(im4, axs[1,1])
	
		im5 = axs[0,2].plot(z_linspace_high, z_dose_high)
		axs[0,2].set_title("z_dose_high")

		im5 = axs[1,2].plot(z_linspace_low, z_dose_low)
		axs[1,2].set_title("z_dose_low")
				
		fig.subplots_adjust(right=0.8)	
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		fig.colorbar(im1, cax=cbar_ax)	
	
		plt.show()


	#print z_linspace		
	return z_dose_low, z_dose_high, z_linspace_low, z_linspace_high

# ########################
# CALCULATE THE RECONSTRUCTED ACTIVITY
# ########################

def calculate_reconstructed_activity ():

	reconstructed_activity = np.fromfile(image_path, dtype=number_format).reshape((mat_x,mat_y,mat_z), order = 'F')
	reconstructed_activity = img_filters.gaussian_filter(reconstructed_activity,(1))
	#print reconstructed_activity.shape
	a=int((fov_x/2-phantom_x/2.)/voxel_x)
	b=int((fov_x/2+phantom_x/2.)/voxel_x)
	c=int((fov_y/2-phantom_y/2.)/voxel_y)
	d=int((fov_y/2+phantom_y/2.)/voxel_y)
	e=int((fov_z/2-phantom_z/2.)/voxel_z)
	f=int((fov_z/2+phantom_z/2.)/voxel_z)
	phantom_reconstructed_activity = reconstructed_activity[ a:b, c:d, e:f ]
	#print phantom_reconstructed_activity.shape
	z_reconstructed_activity = phantom_reconstructed_activity.sum(axis=0).sum(axis=0)
	z_linspace = np.linspace(-phantom_z/2 + voxel_z/2, phantom_z/2-voxel_z/2, activity_dimz)

	plot_flag = 0
	if plot_flag == 1:
		fig, axs = plt.subplots(2,3)

		im1 = axs[0,0].imshow(phantom_reconstructed_activity[:,:,5])
		axs[0,0].set_title("5")
		fig.colorbar(im1, axs[0,0])
		
		im2 = axs[0,1].imshow(phantom_reconstructed_activity[:,:,10])
		axs[0,1].set_title("10")
		fig.colorbar(im2, axs[0,1])
	
		im3 = axs[1,0].imshow(phantom_reconstructed_activity[:,:,15])
		axs[1,0].set_title("15")
		fig.colorbar(im3, axs[1,0])
			
		im4 = axs[1,1].imshow(phantom_reconstructed_activity[:,:,25])
		axs[1,1].set_title("25")
		fig.colorbar(im4, axs[1,1])
	
		im5 = axs[0,2].plot(z_reconstructed_activity)
		axs[0,2].set_title("z_reconstructed_activity")
			
		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		fig.colorbar(im1, cax=cbar_ax)
		
		plt.show()

	#print z_linspace		
		
	return z_reconstructed_activity, z_linspace


# ########################
# CALCULATE THE TRUE ACTIVITY
# ########################

def map_calculate_true_activity(ii):

	beta_map_low_res = np.zeros((activity_dimx, activity_dimy, activity_dimz))
	beta_map_high_res = np.zeros((dose_dimx, dose_dimy, dose_dimz))

	temp_ekine = 0.
	positron_trackID = -1

	gamma_counter = 0
	first_gamma_flag = 0
	temp_eventID = -1

	file_path = proton_beam_path + str(ii) + '.root'
	print file_path
	#try:
	f = ROOT.TFile.Open(file_path)
	#	print 'Analyzed root file number: {0}'.format(i+1)
	#except:
	#	ROOT.TFile.Open(file_path).Close()
	#	print "File {0} doesn't open properly or doesn't exist!".format(file_path)
	#	sys.exit(2)
	
	tree = f.Get("PhaseSpace")
	entries = tree.GetEntries()
	print 'Entries: {0}'.format(entries)
	
	for i, ev in enumerate(tree):
		
		if (i+1)%10000000 == 0:
			print "{0}0 M events analysed".format(int((i+1)/10000000))

		if temp_eventID != ev.EventID:
			temp_eventID = ev.EventID
			
		if (ev.ParticleName.startswith('e+') and (ev.CreatorProcess[0:-1].startswith('RadioactiveDecay'))):
			positron_trackID = ev.TrackID
				
		if (ev.ParticleName.startswith('gamma') and ev.CreatorProcess.startswith('annihil') and ev.ParentID == positron_trackID):
			if not first_gamma_flag:
				temp_ekine = ev.Ekine
				first_gamma_flag = 1

			if first_gamma_flag:

				if (temp_ekine>0.51 and temp_ekine<0.52 and ev.Ekine>0.51 and ev.Ekine<0.512 and \
					ev.EmissionPointX >= -phantom_x/2 and ev.EmissionPointX <= phantom_x/2 and \
					ev.EmissionPointY >= -phantom_y/2 and ev.EmissionPointY <= phantom_y/2 and \
					ev.EmissionPointX >= -phantom_z/2 and ev.EmissionPointZ <= phantom_z/2 ):

					high_x = int((ev.EmissionPointX+phantom_x/2)/dose_voxel_x)
					high_y = int((ev.EmissionPointY+phantom_y/2)/dose_voxel_y)
					high_z = int((ev.EmissionPointZ+phantom_z/2)/dose_voxel_z)
					if (ev.EmissionPointX == phantom_x/2.):
						print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
						high_x -= 1
					if (ev.EmissionPointY == phantom_y/2.):
						print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
						high_y -= 1
					if (ev.EmissionPointZ == phantom_z/2.):
						print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
						high_z -= 1

					low_x = int((ev.EmissionPointX+phantom_x/2)/voxel_x)
					low_y = int((ev.EmissionPointY+phantom_y/2)/voxel_y)
					low_z = int((ev.EmissionPointZ+phantom_z/2)/voxel_z)
					if (ev.EmissionPointX == phantom_x/2.):
						print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
						low_x -= 1
					if (ev.EmissionPointY == phantom_y/2.):
						print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
						low_y -= 1
					if (ev.EmissionPointZ == phantom_z/2.):
						print "Source position at the boundary: {0}, {1}, {2}".format(entry.sourcePosX1, entry.sourcePosY1, entry.sourcePosZ1)
						low_z -= 1

					beta_map_low_res [low_x, low_y, low_z] +=1
					beta_map_high_res [high_x, high_y, high_z] +=1

					gamma_counter += 1
				first_gamma_flag = 0
				
	np.save(positron_map_low_path+str(ii), beta_map_low_res)
	np.save(positron_map_high_path+str(ii), beta_map_high_res)

	print 'Gammas counter for {0} phase space file: {1}'.format(ii, gamma_counter)



def calculate_true_activity ():

	beta_map_low_res = np.zeros((activity_dimx, activity_dimy, activity_dimz))
	beta_map_high_res = np.zeros((dose_dimx, dose_dimy, dose_dimz))

	no_root_files = 20;

	maps_list = range(1, no_root_files+1)
	print maps_list
	p = Pool(20)
	#p.map (map_calculate_true_activity, maps_list)
	#map_calculate_true_activity(1)


	for i in range(1, no_root_files+1):
		print 'Load beta maps from file: {0}'.format(i)
		temp_high=np.load(positron_map_high_path+str(i)+'.npy')
		temp_low=np.load(positron_map_low_path+str(i)+'.npy')
		beta_map_high_res = beta_map_high_res + temp_high
		beta_map_low_res = beta_map_low_res + temp_low


	true_z_low = beta_map_low_res.sum(axis=0).sum(axis=0)
	true_z_high = beta_map_high_res.sum(axis=0).sum(axis=0)

	true_linspace_low = np.linspace(-phantom_z/2 + voxel_z/2, phantom_z/2-voxel_z/2, activity_dimz)		
	true_linspace_high = np.linspace(-phantom_z/2 + dose_voxel_z/2, phantom_z/2-dose_voxel_z/2, dose_dimz)		

	plot_flag = 0
	if plot_flag == 1:

		fig, axs = plt.subplots(2,3)

		im1 = axs[0,0].imshow(beta_map_high_res[:,:,24])
		axs[0,0].set_title("24")
		fig.colorbar(im1, axs[0,0])
			
		im2 = axs[0,1].imshow(beta_map_high_res[:,:,73])
		axs[0,1].set_title("73")
		fig.colorbar(im2, axs[0,1])
	
		im3 = axs[1,0].imshow(beta_map_high_res[:,:,74])
		axs[1,0].set_title("74")
		fig.colorbar(im3, axs[1,0])
			
		im4 = axs[1,1].imshow(beta_map_high_res[:,:,75])
		axs[1,1].set_title("75")
		fig.colorbar(im4, axs[1,1])
	
		im5 = axs[0,2].plot(true_linspace_high, true_z_high)
		axs[0,2].set_title("true_z_high")

		im5 = axs[1,2].plot(true_linspace_low, true_z_low)
		axs[1,2].set_title("true_z_low")
				
		fig.subplots_adjust(right=0.8)	
		cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
		fig.colorbar(im1, cax=cbar_ax)	
	
		plt.show()

	return true_z_low, true_z_high, true_linspace_low, true_linspace_high
	
# ########################
# MAIN
# ########################


if __name__ == '__main__':



	print 'Calculating dose ...'
	dose_low, dose_high, linspace_dose_low, linspace_dose_high = calculate_dose()
	print 'Calculating reconstructed beta profile ...'
	reconstructed, linspace_recon = calculate_reconstructed_activity()
	print 'Calculating true beta profile ...'
	true_low, true_high, linspace_true_low, linspace_true_high = calculate_true_activity()
	#calculate_true_activity()
		

	dose = dose_high
	x_dose = linspace_dose_high
	true = true_high
	x_true = linspace_true_high
	
	low = 1	
	if low:
		dose = dose_low
		x_dose = linspace_dose_low
		true = true_low
		x_true = linspace_true_low


	fig, ax1 = plt.subplots()

	ax1.plot(linspace_recon, reconstructed, color='b', linestyle = '--', linewidth = 3.)
	ax1.plot(x_true, true, color='g', linestyle = ':', linewidth = 3.)
	ax1.set_ylim([1,1000000])
	if low:
		ax1.set_ylim([1,10000000])
	ax1.set_xlabel('Phantom position [mm]', fontweight="bold", fontsize = 20)
	ax1.set_ylabel('Activity [counts]', fontweight="bold", fontsize = 20)
	ax1.set_title('Reconstructed activity and deposited dose for the proton beam irradiation\nof a 5x5x20cm3 PMMA phantom', fontweight="bold", fontsize = 24)
	ax1.grid(linestyle = ':', linewidth = 0.75, color = 'grey')
 	ax1.tick_params('x', width=2, labelsize=20)   
	ax1.tick_params('y', width=2, labelsize=20)   
	ax1.set_yscale('log')

	ax2 = ax1.twinx()    
	ax2.plot(x_dose, dose, color = 'magenta', linestyle = '-', linewidth = 3.)
	ax2.set_ylabel('Dose [A.U.]', fontweight="bold", color='magenta', fontsize = 20)
	ax2.tick_params('y', colors='magenta', width=2, labelsize=20) 
	ax1.legend(['Reconstructed profile', 'MC activity profile'], loc=2, fontsize = 'x-large')
	ax2.legend(['Deposited dose'], fontsize = 'x-large')

	plt.show()

#	plt.plot (iterr, nrmsd, 'ro')
#	plt.axis([0, 100, 0.975, 0.98])

#	plt.title('NRMSD calculated based on cubic water phantom simulation with 1G primary back-to-back source uniformly distributed within the phantom \n(Gaitanis, Anastasios, et al. "PET image reconstruction: A stopping rule for the MLEM algorithm based on properties of the updating coefficients."\n Computerized Medical Imaging and Graphics 34.2 (2010): 131-141.)' , fontweight="bold")

#	plt.xlabel('Iteration number' , fontweight="bold")
#	plt.ylabel('NRMSD' , fontweight="bold")

#	plt.grid(True)

#	plt.show()



