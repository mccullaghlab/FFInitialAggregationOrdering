import os
import numpy as np 
from sklearn.decomposition import PCA
from scipy.linalg import logm, expm
import MDAnalysis
import matplotlib.pyplot as plt





def compute_pbc_dist2(r1,r2,box,hbox):
        dist2 = 0
        for j in range(0,3):
            temp = r1[j]-r2[j]
            if temp < -hbox[j]:
                temp += box[j]
            elif temp > hbox[j]:
                temp -= box[j]
            dist2 += temp*temp
        return dist2;
def computeOrthogVec(r1,r2,r3,box): 
        vec1=np.zeros(3, dtype=np.float)
        vec2=np.zeros(3, dtype=np.float)
        for i in range(0,3):
            vec1[i] = r2[i] - r1[i]
            vec2[i] = r3[i] - r2[i]
        crossvec = np.cross(vec1,vec2)
        crossvec /= np.linalg.norm(crossvec)
        return crossvec; 

traj_file = "RemoveMonomerFromCrystal.w??.runXX_nowater.xyz.nc"
top_file = "RemoveMonomerFromCrystal.stripwater.prmtop"
coord = MDAnalysis.Universe(top_file, traj_file)

# make an atom selection
NT=['resid 15 and name N H1 H2 H3']
CT=['resid 16 and name C O OXT']

count=-1
coord.atoms.n_residues
data= np.zeros((len(coord.trajectory),3))
for ts in coord.trajectory[0::10]:
   count+=1
   box = ts.dimensions[0:3]
   hbox = box/2.0
   NT= coord.select_atoms('resid 15 and name N H1 H2 H3').center_of_mass()
   CT= coord.select_atoms('resid 16 and name C O OXT').center_of_mass()
   mon1=coord.select_atoms('resid 15 16 and backbone').center_of_mass()
   mon2=coord.select_atoms('resid 83 and name C or resid 84 and name N').center_of_mass()
   dist2=compute_pbc_dist2(mon1,mon2,box,hbox)
   dist=np.sqrt(dist2)
   vector = NT-CT
   zvector = vector[2]/np.linalg.norm(vector) 

   data[(count,0)] = dist
   data[(count,1)] = zvector
   data[(count,2)] = vector[2]
   #data[(count,1)] = np.arccos(zvector)*180./np.pi


np.savetxt("AlignmentOfMonomerw??_runXX.dat",data)
