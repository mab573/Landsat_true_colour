#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# SCRIPT NAME		:: Make_final_lookups.py 

# PURPOSE 		:: Create the final lookup files that will be used to perform atmospheric correction    

# SYNOPSIS 		:: Open a bunch of MODTRAN .7sc files and calculate parameters     

# MODULES CALLED 	:: numpy
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

import matplotlib
from numpy import *
import matplotlib.pyplot as plt
import glob
import os,sys,getopt

'''
#*******************************************************

def usage():

 print "SYNOPSIS:\n This program creates lookup files that are used in the atmospheric correction process"
 print "\nUSAGE: Make_final_lookups.py [OPTIONS] ";
 print "\nOPTIONS:\n";
 print "\n-i, (Working directory)";
 print "Base Usage:"
 print "python Make_final_lookups.py -i /Dir/sub_dir\n\n"
 sys.exit()

#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>********************
#** Read in command line arguments like input hyperion L1R file, auxillary file and header file 
#** Use getopt so that all input arguments are strictly defined by switch and not position in argument list

try:
	opts, args = getopt.getopt(sys.argv[1:], "hi:", ["help", "in_dir"])

except getopt.GetoptError, err:

        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

for o, a in opts:
	if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--in_dir"):
            in_dir = a
	else:
            assert False, "unhandled option"
    

#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>********************

'''

def Make_final_lookups(in_dir):

	##############################################################################################################################################################
	#
	# The rho1_array and rho2_array are read in columns from MODTRAN5 tape 7 files. 
	# The two files are supposed to be exactly the same atmospheric conditions, observation geometry and so on but for two different surface reflectance values.
	# Using the process outlined in Guanter et al. and the little document from Andrew, this information an be used to perform atmospheric compensation/correction.

	#These columns are:
	#FREQ  TOT_TRANS  PTH_THRML  THRML_SCT  SURF_EMIS   SOL_SCAT  SING_SCAT  GRND_RFLT  DRCT_RFLT  TOTAL_RAD  REF_SOL  SOL@OBS   DEPTH DIR_EM    TOA_SUN BBODY_T[K]

	## Constants

	rho1=0.50	#may want to have these as inputs
	rho2=1.00

	## READ IN THE TAPE7 FILES ##
	## IDENTIFY THE PAIRS OF 50% AND 100% FILES

	files_50_percent=[]
	files_100_percent=[]

	os.chdir(in_dir) 
	for files in glob.glob("*_0.50.7sc"):
		files_50_percent.append(files)

	for files in glob.glob("*_1.00.7sc"):
		files_100_percent.append(files)

	#print files_50_percent
	#print ' ' 
	#print files_100_percent
	#sys.exit()

	####
	####
	####
	#### LOOP THROUGH FILES AND CREATE A SET OF LOOKUP FILES FOR A RANGE OF WATER VAPOUR AMOUNTS

	for xx in range(len(files_50_percent)):

		rho1_tape7_array=loadtxt(files_50_percent[xx], dtype=float, skiprows = 11, comments='-9999.')
		rho2_tape7_array=loadtxt(files_100_percent[xx], dtype=float, skiprows = 11, comments='-9999.')


		## Assign tape7 columns to the correct terms

		t_dir_up_rho1 = rho1_tape7_array[:,1]	#this is exactly the same as t_dir_up_rho2
		L_ss_rho1 = rho1_tape7_array[:,5]
		L_gr_rho1 = rho1_tape7_array[:,7]
		L_tr_rho1 = rho1_tape7_array[:,9]

		t_dir_up_rho2 = rho2_tape7_array[:,1]
		L_ss_rho2 = rho2_tape7_array[:,5]
		L_gr_rho2 = rho2_tape7_array[:,7]
		L_tr_rho2 = rho2_tape7_array[:,9]

		## Calculate wavelength in nanometers

		wavelength=rho1_tape7_array[:,0]

		## Calculate downwelling flux (direct and diffuse)

		Eg_rho1=(math.pi*L_gr_rho1)/(t_dir_up_rho1*rho1)
		Eg_rho2=(math.pi*L_gr_rho2)/(t_dir_up_rho2*rho2)

		## Calculate the path radiance

		Lp_0=((rho1*Eg_rho1*L_ss_rho2)-(rho2*Eg_rho2*L_ss_rho1))/((rho1*Eg_rho1)-(rho2*Eg_rho2))

		## Calculate the diffuse transmittance

		t_dif=(math.pi*(L_ss_rho1-Lp_0))/(rho1*Eg_rho1)

		## Calculate the spherical albedo

		S=(Eg_rho2-Eg_rho1)/((rho2*Eg_rho2)-(rho1*Eg_rho1))

		## Calculate the direct and diffuse upward transmittance

		Gamma_up=t_dif+t_dir_up_rho1

		## Calculate downwelling flux (direct and diffuse) decoupled from the surface

		Eg_0=(math.pi*(L_tr_rho1-Lp_0)*(1-(S*rho1)))/(rho1*Gamma_up)

		## Calculate cosine of the solar zenith angle

		# mu_s_2=(math.pi*rho2_tape7_array[:,8])/(rho2*(rho2_tape7_array[:,10])) # These aren't used for anything so no need to calculate them
		# mu_s_1=(math.pi*rho1_tape7_array[:,8])/(rho1*(rho1_tape7_array[:,10]))


		## Replace nan values (caused by divide by zero in above calcs) with a number mainly for text formating

		for n,i in enumerate(Eg_0):
			if math.isnan(i) == 1:
				Eg_0[n]=-9.9999e-99

		for n,i in enumerate(Lp_0):
			if math.isnan(i) == 1:
				Lp_0[n]=-9.9999e-99

		for n,i in enumerate(Gamma_up):
			if math.isnan(i) == 1:
				Gamma_up[n]=-9.9999e-99

		for n,i in enumerate(S):
			if math.isnan(i) == 1:
				S[n]=-9.9999e-99

		## Write required calculated values to a file as part of a lookup table
	
		out_50=files_50_percent[xx]
		out_base_name_50=out_50.split('*.7sc')

		savetxt(out_base_name_50[0]+'_lookup.txt', transpose([Eg_0,Lp_0,Gamma_up,S]),fmt="%1.4e", delimiter='\t')

		#print "xx =", xx, out_50

	## This will be the same for all files so it only needs to be saved once
	savetxt('lookup_wavelengths.txt', wavelength,fmt="%4.2f") 
