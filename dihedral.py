#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# USAGE:
# ./rgyr.py system_descripter pdb_file traj_file 


# PREAMBLE:

import numpy as np
import MDAnalysis
import sys
import math

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

def computeDihedral(r1, r2, r3, r4, box):
#####Praxeolitic formula#####
	vec1=np.zeros(3, dtype=np.float)
	vec2=np.zeros(3, dtype=np.float)
	vec3=np.zeros(3, dtype=np.float)
 
	for j in range(0,3):
		# Define the 3 vectors within the dihedral
		vec1[j] = -1.0*(r2[j]-r1[j])
		vec2[j] = r3[j]-r2[j]
		vec3[j] = r4[j]-r3[j]

		# Pbc
		#if vec1[j] < -box[j]/2.0:
                #        vec1[j] += box[j]
		#elif vec1[j] > box[j]/2.0:
                #        vec1[j] -= box[j]
		#if vec2[j] < -box[j]/2.0:
                #        vec2[j] += box[j]
		#elif vec2[j] > box[j]/2.0:
                #        vec2[j] -= box[j]
		#if vec3[j] < -box[j]/2.0:
                #        vec3[j] += box[j]
		#elif vec3[j] > box[j]/2.0:
                #        vec3[j] -= box[j]
        
        vec2 /= np.linalg.norm(vec2)
	v = vec1 - np.dot(vec1, vec2)*vec2
        w = vec3 - np.dot(vec3, vec2)*vec2
        x = np.dot(v, w)
        y = np.dot(np.cross(vec2, v), w)
        dih = np.arctan2(y, x)
        ## Calculate angle between 2 planes
	#n1_v = np.cross(vec1, vec2)
	#n2_v = np.cross(vec2, vec3)
	#m_v = np.cross(n1_v, n2_v)
	##m_v = np.cross(n1_v, vec2)
	#
	#n1_norm = n1_v/np.linalg.norm(n1_v)
	#n2_norm = n2_v/np.linalg.norm(n2_v)
	#m_norm = m_v/np.linalg.norm(m_v)
	#
	#x = np.dot(n1_norm, n2_norm)
	#y = np.dot(m_norm, n2_norm)
	#dih = -np.arctan2(y,x) 

	return dih;

system = sys.argv[1]
pdb = sys.argv[2]
traj = sys.argv[3]

#res_list = [2]
res_list = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120]
#res_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120]
omega_sel = []
phi_sel = []
psi_sel = []
jake_sel = []
chi1I_sel = []
chi2I_sel = []
chi1II_sel = []
chi2II_sel = []

u = MDAnalysis.Universe('%s' %(pdb), '%s' %(traj))

for k in range(len(res_list)):
#	#for even resids
#	if res_list[k] % 2 == 0:
	i = res_list[k] - 1
	j = res_list[k]
#	#for odd resids
#	if res_list[k+1] % 2 == 0:
#		i = res_list[k]
#		j = i + 1
#		l = i - 1		
		
#	#else:
	omega_sel1 = []
	omega_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
	omega_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
	omega_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
	omega_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
	omega_sel.append(omega_sel1)
	
	phi_sel1 = []
	phi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
	phi_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
	phi_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
	phi_sel1.append(u.select_atoms("resid %s and name C" %(j)))  
	phi_sel.append(phi_sel1)
	
	psi_sel1 = []
        psi_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
        psi_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
        psi_sel1.append(u.select_atoms("resid %s and name C" %(i)))  
        psi_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
        psi_sel.append(psi_sel1)

	jake_sel1 = []
	jake_sel1.append(u.select_atoms("resid %s and name CZ" %(i)))  
	jake_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
	jake_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
	jake_sel1.append(u.select_atoms("resid %s and name CZ" %(j)))  
	jake_sel.append(jake_sel1)
	
	chi1I_sel1 = []
	chi1I_sel1.append(u.select_atoms("resid %s and name N" %(i)))  
	chi1I_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
	chi1I_sel1.append(u.select_atoms("resid %s and name CB" %(i)))  
	chi1I_sel1.append(u.select_atoms("resid %s and name CG" %(i)))  
	chi1I_sel.append(chi1I_sel1)

	chi2I_sel1 = []
	chi2I_sel1.append(u.select_atoms("resid %s and name N" %(j)))  
	chi2I_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
	chi2I_sel1.append(u.select_atoms("resid %s and name CB" %(j)))  
	chi2I_sel1.append(u.select_atoms("resid %s and name CG" %(j)))  
	chi2I_sel.append(chi2I_sel1)
	
	chi1II_sel1 = []
	chi1II_sel1.append(u.select_atoms("resid %s and name CA" %(i)))  
	chi1II_sel1.append(u.select_atoms("resid %s and name CB" %(i)))  
	chi1II_sel1.append(u.select_atoms("resid %s and name CG" %(i)))  
	chi1II_sel1.append(u.select_atoms("resid %s and name CD1" %(i)))  
	chi1II_sel.append(chi1II_sel1)

	chi2II_sel1 = []
	chi2II_sel1.append(u.select_atoms("resid %s and name CA" %(j)))  
	chi2II_sel1.append(u.select_atoms("resid %s and name CB" %(j)))  
	chi2II_sel1.append(u.select_atoms("resid %s and name CG" %(j)))  
	chi2II_sel1.append(u.select_atoms("resid %s and name CD1" %(j)))  
	chi2II_sel.append(chi2II_sel1)
	
	#

	#print k,phi_sel[k][0],phi_sel[k][1],phi_sel[k][2],phi_sel[k][3]
	#print k,psi_sel[k][0],psi_sel[k][1],psi_sel[k][2],psi_sel[k][3]
	#print k,omega_sel[k][0],omega_sel[k][1],omega_sel[k][2],omega_sel[k][3]
#
#	if res_list[k] =

#	phi_sel.append(u.select_atoms("resid %s and name C" %(i)), u.select_atoms("resid %s and name CA" %(i)), u.select_atoms("resid %s and name N" %(i)), u.select_atoms("resid %s and name C" %(j)))
#	psi_sel.append(u.select_atoms("resid %s and name N" %(i)), u.select_atoms("resid %s and name C" %(j)), u.select_atoms("resid %s and name CA" %(j)), u.select_atoms("resid %s and name N" %(j)))
	
#sel1=u.select_atoms('resid 1:13 and resname PDI')
#sel2=u.select_atoms('resid 14:26 and resname PDI')
#sel3=u.select_atoms('resid 11').residues[0].phi_selection()
#sel4=u.select_atoms('resid 11').residues[0].psi_selection()
#sel5=u.select_atoms('resid 11').residues[0].omega_selection()
#print sel3, sel4, sel5
# SUBROUTINES:

# MAIN PROGRAM:

out1 = open('%s.omega-dihedral.dat' %(system), 'w')
out2 = open('%s.phi-dihedral.dat' %(system), 'w')
out3 = open('%s.psi-dihedral.dat' %(system), 'w')
out4 = open('%s.CZ-CA-CA-CZdihedral.dat' %(system), 'w')
out5 = open('%s.oddChi1.dat' %(system), 'w')
out6 = open('%s.evenChi1.dat' %(system), 'w')
out7 = open('%s.oddChi2.dat' %(system), 'w')
out8 = open('%s.evenChi2.dat' %(system), 'w')

omega_data = np.zeros(len(res_list), dtype=np.float)
phi_data = np.zeros(len(res_list), dtype=np.float)
psi_data = np.zeros(len(res_list), dtype=np.float)
jake_data = np.zeros(len(res_list), dtype=np.float)
chi1I_data = np.zeros(len(res_list), dtype=np.float)
chi2I_data = np.zeros(len(res_list), dtype=np.float)
chi1II_data = np.zeros(len(res_list), dtype=np.float)
chi2II_data = np.zeros(len(res_list), dtype=np.float)

for ts in u.trajectory:
	for i in range(len(res_list)):
	        #	print i,phi_sel[i][0],phi_sel[i][1],phi_sel[i][2],phi_sel[i][3]
	       	omega_data[i]=computeDihedral(omega_sel[i][0].center_of_mass(), omega_sel[i][1].center_of_mass(), omega_sel[i][2].center_of_mass(), omega_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	phi_data[i]=computeDihedral(phi_sel[i][0].center_of_mass(), phi_sel[i][1].center_of_mass(), phi_sel[i][2].center_of_mass(), phi_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	psi_data[i]=computeDihedral(psi_sel[i][0].center_of_mass(), psi_sel[i][1].center_of_mass(), psi_sel[i][2].center_of_mass(), psi_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	jake_data[i]=computeDihedral(jake_sel[i][0].center_of_mass(), jake_sel[i][1].center_of_mass(), jake_sel[i][2].center_of_mass(), jake_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	chi1I_data[i]=computeDihedral(chi1I_sel[i][0].center_of_mass(), chi1I_sel[i][1].center_of_mass(), chi1I_sel[i][2].center_of_mass(), chi1I_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	chi2I_data[i]=computeDihedral(chi2I_sel[i][0].center_of_mass(), chi2I_sel[i][1].center_of_mass(), chi2I_sel[i][2].center_of_mass(), chi2I_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	chi1II_data[i]=computeDihedral(chi1II_sel[i][0].center_of_mass(), chi1II_sel[i][1].center_of_mass(), chi1II_sel[i][2].center_of_mass(), chi1II_sel[i][3].center_of_mass(), u.dimensions[:3])
	       	chi2II_data[i]=computeDihedral(chi2II_sel[i][0].center_of_mass(), chi2II_sel[i][1].center_of_mass(), chi2II_sel[i][2].center_of_mass(), chi2II_sel[i][3].center_of_mass(), u.dimensions[:3])


	        out1.write(' %10.6f' %(omega_data[i]*180/np.pi))
	        out2.write(' %10.6f' %(phi_data[i]*180/np.pi))
	        out3.write(' %10.6f' %(psi_data[i]*180/np.pi))
	        out4.write(' %10.6f' %(jake_data[i]*180/np.pi))
	        out5.write(' %10.6f' %(chi1I_data[i]*180/np.pi))
	        out6.write(' %10.6f' %(chi2I_data[i]*180/np.pi))
	        out7.write(' %10.6f' %(chi1II_data[i]*180/np.pi))
	        out8.write(' %10.6f' %(chi2II_data[i]*180/np.pi))
		
#		if i == 5:
#			dih_phi = sel3.dihedral.value()
#			dih_psi = sel4.dihedral.value()
#			dih_omega = sel5.dihedral.value()
#			print dih_phi, phi_data[i]*180/np.pi
#			print dih_psi, psi_data[i]*180/np.pi
#			print dih_omega, omega_data[i]*180/np.pi
			
		
	out1.write('\n')
	out2.write('\n')
	out3.write('\n')
	out4.write('\n')
	out5.write('\n')
	out6.write('\n')
	out7.write('\n')
	out8.write('\n')
out1.close()
out2.close()
out3.close()
out4.close()
out5.close()
out6.close()
out7.close()
out8.close()

