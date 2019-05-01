#!/usr/bin/env python

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# SCRIPT NAME		:: Create_Tape5_files.py 

# PURPOSE 		:: This script will take a template MODTRAN5 tape 5 file as input and generate a series of tape 5 files with a range of
#			:: water vapour values. This tape 5 file will have to come from whatever creates the tape 5 files for the lookups
#			:: which are used to do the atmospheric correction. All of the geometry and other parameters must be the same.    

# SYNOPSIS 		:: The template file will be read in. A bunch of tape 5 files will be created and run with a range of wv values.
#			:: These will extract band averaged values for three band, one in a wv absoption feature and on either side.
#			:: These will be used to calculate a water vapour value based on a band ratio.     

# MODULES CALLED 	:: numpy
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#**************** Load modules and such ***************

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import getopt, sys, math, Image, os, glob, subprocess
import csv

#*******************************************************

'''

def usage():

 print "SYNOPSIS:\n "
 print "\nUSAGE:  ";
 print "\nOPTIONS:\n";
 print "\n-i, Day of year extrcted from the metadata file";
 print "\n\n"
 print "\n-m, Scene centre lat and lon estimated from the corner lat lons from the metadata file.";
 print "\n\n"
 print "\n-e, Time of the overpass extracted from the metadata file";
 print "\n\n"
 print "\n-t, The path to the tape5 directory.";
 print "\n\n"
 print "\n-t, The path to the MODTRAN data directory.";
 print "\n-h, print the usage statement"; 
 print "Base Usage:"
 print "\n\n"
 sys.exit()


##############################################################################################################################################################
#
# 

## Constants

#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>********************
#** Read in command line arguments like input hyperion L1R file, auxillary file and header file 
#** Use getopt so that all input arguments are strictly defined by switch and not position in argument list

try:
	opts, args = getopt.getopt(sys.argv[1:], "hi:m:e:t:d:", ["help", "day_of_year", "scene_centre","op_time","tape5_dir","data_dir"])

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
        elif o in ("-i", "--day_of_year"):
            day_of_year = a
	elif o in ("-m", "--scene_centre"):
            scene_centre = a
	elif o in ("-e", "--op_time"):
            op_time = a
	elif o in ("-t", "--tape5_dir"):
            tape5_dir = a
	elif o in ("-d", "--data_dir"):
            data_dir = a
	else:
            assert False, "unhandled option"
    

#***************<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>********************

'''

#********************************************
## This reads the number of lines in a file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
#*******************************************


def write_tape5(day_of_year,scene_centre, op_time, tape5_dir, data_dir):


	#*******************************************************************************
	# LIST THE VARIABLE PARAMETERS TO BE APPLIED TO THE GENERATED TAPE5 FILES ******

	# After some investigation it seems that small increases of wv have a large
	# effect on the depth of the feature at the lower end of the scale. Its probably
	# better to have more resolution at this end.

	'''
	CWV=['0.05','0.08','0.10','0.12','0.14','0.18','0.20','0.25',
	'0.30','0.35','0.40','0.45','0.50','0.60','0.70','0.80',
	'0.90','1.00','1.10','1.20','1.30','1.40','1.50','1.60',
	'1.70','1.85','2.10','2.35','2.60','3.00','3.50','4.00',
	'4.50','5.00']
	'''

	look_angle='180.000'

	CWV=['0.10']

	REF=['0.50','1.00']

	in_tape5_file='/Users/212062O/src/python/Landsat_Processing/template_files/MOD_Template_Landsat.tp5'

	#*******************************************************************************

	## READ IN THE BASE TAPE5 FILE AND SAVE IT TO THE LIST LINES ##
	lines=[]
	new_file_name=[]



	if in_tape5_file is None:

    		sys.exit()

	else:

		with open(in_tape5_file,'rb') as f:

			for line in f.readlines():
		
				lines.append(line)
		
		f.closed


	## CHANGE ALL OF THE STUFF IN THE TAPE5 FILE FOR EACH SCENE THAT IS NOT INCREMENTED

	new_line=lines[4]
	new_line=new_line[:23]+look_angle+new_line[30:]
	lines[4]=new_line

	new_line=lines[5]
	new_line=new_line[:12]+day_of_year+new_line[15:]
	lines[5]=new_line

	new_line=lines[6]
	new_line=new_line[:3]+str(scene_centre[0])+new_line[10:]
	new_line=new_line[:12]+str(scene_centre[1])+new_line[19:]
	new_line=new_line[:44]+op_time+new_line[50:]
	lines[6]=new_line

	# Open the modroot.in file so that all of the tape5 file names can be written to it.



	modroot = open('mod5root.in', "w")


	## SET UP LOOPS TO INDEX THE LISTS/ARRAYS WITH THE VARIABLE MODTRAN PARAMETERS ##

	CWV_len=len(CWV)
	file_lines=file_len(in_tape5_file)

	print lines[1]

	new_ref_line=lines[0]
	new_wv_line=lines[1]


	for REF_cnt in range(0,2):

		new_ref_line=new_ref_line[:74]+REF[REF_cnt]+new_ref_line[78:]
		lines[0]=new_ref_line
		print REF_cnt

		for CWV_cnt in range(0,CWV_len):

			# The WV amount is set in card 1A which is required. It should always be the second line in the tape 5 file and it should always be in the same spot.
			new_wv_line=new_wv_line[:23]+CWV[CWV_cnt]+new_wv_line[27:]
			lines[1]=new_wv_line
	
 			new_file_name=tape5_dir+'MOD_ATM-cm_WV_'+CWV[CWV_cnt]+'_'+REF[REF_cnt]+'.tp5'

			# Write the new file name to the modroot.in file including the full path

			modroot_entry=tape5_dir+'MOD_ATM-cm_WV_'+CWV[CWV_cnt]+'_'+REF[REF_cnt]+'.tp5\n'
			modroot.writelines(modroot_entry)

			# Open new file to write stuff to
			new_file = open(new_file_name, "w")

			for xx in range(file_lines):
	
				new_file.writelines(lines[xx])

			new_file.close()

	modroot.close()

	## To run MODTRAN from any folder it needs to contain the mod5root.in file and a symbolic link called DATA
	## pointing to the DATA directory from the MODTRAN installation. The following command creates the symbolic link.

	subprocess.call(["ln", "-s", data_dir,"DATA"])

