#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
#USAGE :  python python_mdanalysis_skeleton.py [config file name]                                                                  

#DEPENDCIES : numpy, MDAnalysis, math, sys                                                                                         

#CONFIG FILE FORMAT:                                                                                                               
#   TopFile = [topology file name (prmtop file)]                                                                                   
#   TrajFile = [trajectory file name (mdcrd file)]                                                                                 
#   OutFile = [output data file name]                                                                                              

import numpy as np
import sys
import os
import MDAnalysis
import math
# compute the distance between two points taking into account periodic boundary conditions                                         
def computePbcDist2(r1,r2,box):
        dist2 = 0

        for j in range(0,3):
                temp = r1[j]-r2[j]
                if temp < -box[j]/2.0:
                        temp += box[j]
                elif temp > box[j]/2.0:
                        temp -= box[j]
                dist2 += temp*temp

        return dist2;
# read the configuration file and populate the global variables                                                                    
def ParseConfigFile(cfg_file):
        global top_file, traj_file, out_file
        f = open(cfg_file)
        for line in f:
                # first remove comments                                                                                            
                if '#' in line:
                        line, comment = line.split('#',1)
                if '=' in line:
                        option, value = line.split('=',1)
                        option = option.strip()
                        value = value.strip()
                        print "Option:", option, " Value:", value
                        # check value                                                                                              
                        if option.lower()=='topfile':
                                top_file = value
                        elif option.lower()=='trajfile':
                                traj_file = value
                        elif option.lower()=='outfile':
                                out_file = value
                        else :
                                print "Option:", option, " is not recognized"
        f.close()

# Main Program                                                                                                                     

# read in command line argument                                                                                                    
cfg_file = sys.argv[1]

# read cfg file                                                                                                                    
ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Trajectory file:", traj_file
print "Output data file:", out_file

# initiate MDAnalysis coordinate universe                                                                                          
coord = MDAnalysis.Universe(top_file, traj_file)
#coord = MDAnalysis.Universe(traj_file, format="ncdf")
# make an atom selection                                                                                                           
sel_list1=[]
sel_list2=[]
odd_residues=[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65,  67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119] 
nCG_sites=60
for i in odd_residues:
        j=i+1
        temp1=coord.select_atoms('resid %s:%s and name C CA N O OXT HA H H1 H3 H2'%(i,j)) #backbone
        #temp1=coord.select_atoms('resid %s:%s and name CB HB2 HB3 CG CD1 HD1 CE1 HE1 CZ HZ CE2 HE2 CD2 HD2'%(i,j)) #sidechain
        #temp1=coord.select_atoms('resid %s and name CB HB2 HB3 CG CD1 HD1 CE1 HE1 CZ HZ CE2 HE2 CD2 HD2'%(i)) #odd phenyl ring
        #temp2=coord.select_atoms('resid %s and name CB HB2 HB3 CG CD1 HD1 CE1 HE1 CZ HZ CE2 HE2 CD2 HD2'%(j)) #even phenyl ring
	sel_list1.append(temp1)
	#sel_list2.append(temp2)
# declare arrays of histogram                                                                                                      
hist_min = 0.
hist_max = 30.
bin_size = 0.1
num_bins = int((hist_max-hist_min)/bin_size)
dist_hist = np.zeros(num_bins,dtype=float)
box = coord.dimensions[:3]
count = 0
print len(sel_list1)
volume = box[0]*box[1]*box[2]
# Loop through trajectory                                                                                                          
for ts in coord.trajectory:
	for i in range(0, nCG_sites -1):
		for j in range(i+1, nCG_sites):
                        count+=1
			dist1 = math.sqrt(computePbcDist2(sel_list1[i].center_of_mass(),sel_list1[j].center_of_mass(),coord.dimensions[:3]))
			#dist2 = math.sqrt(computePbcDist2(sel_list1[i].center_of_mass(),sel_list2[j].center_of_mass(),coord.dimensions[:3]))
			#dist3 = math.sqrt(computePbcDist2(sel_list2[i].center_of_mass(),sel_list2[j].center_of_mass(),coord.dimensions[:3]))
			dist_bin1 = int((dist1-hist_min)/bin_size)
			#dist_bin2 = int((dist2-hist_min)/bin_size)
			#dist_bin3 = int((dist3-hist_min)/bin_size)
			if hist_min<dist1<hist_max:
				dist_hist[dist_bin1] += 1
			#if hist_min<dist2<hist_max:
			#	dist_hist[dist_bin2] += 1
			#if hist_min<dist3<hist_max:
			#	dist_hist[dist_bin3] += 1

#Normalize                                                                      
for i in range(num_bins):
        dist_hist[i] /= (count*4.*math.pi*bin_size*((i+0.5)*bin_size+hist_min)**2)/volume

# open output files                                                                                                                
out = open(out_file,'w')

for i in range(num_bins):
        out.write("%10.5f %10.5f\n" %(i*bin_size+hist_min,dist_hist[i]))

# close output file                                                                                                                
out.close
