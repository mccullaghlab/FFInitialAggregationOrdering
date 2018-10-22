#USAGE :  python [python script] [config file name]                                                                  


#CONFIG FILE FORMAT:                                                                                                               
#   TopFile = [topology file name (prmtop or pdb(if pdb makesure it has the same ordering of atoms as the trajectories topology file. If it does not, any results will be incorrect))]                      
#   TrajFile = [trajectory file name (.nc)]                                            
#   OutFiles = [output data file name]                                                                                              

import numpy as np
import sys
import os
import MDAnalysis
#from MDAnalysis.tests.datafiles import DCD
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
        global top_file, traj_file, out1_file, outtraj1, outtraj2, outtraj3, outtraj4, outtraj5, outtraj6, outtraj7, outtraj8, outtraj9, outtraj10, outtraj11, outtraj12, outtraj13, outtraj14
        f = open(cfg_file)
        for line in f:
                # first remove comments                                                                                            
                if '#' in line:
                        line, comment = line.split('#',1)
                if '=' in line:
                        option, value = line.split('=',1)
                        option = option.strip()
                        value = value.strip()
                        #print "Option:", option, " Value:", value
                        # check value                                                                                              
                        if option.lower()=='topfile':
                                top_file = value
                        elif option.lower()=='trajfile':
                                traj_file = value
                        elif option.lower()=='outfile1':
                                out1_file = value
                        elif option.lower()=='outfile2':
                                outtraj1 = value
                        elif option.lower()=='outfile3':
                                outtraj2 = value
                        elif option.lower()=='outfile4':
                                outtraj3 = value
                        elif option.lower()=='outfile5':
                                outtraj4 = value
                        elif option.lower()=='outfile6':
                                outtraj5 = value
                        elif option.lower()=='outfile7':
                                outtraj6 = value
                        elif option.lower()=='outfile8':
                                outtraj7 = value
                        elif option.lower()=='outfile9':
                                outtraj8 = value
                        elif option.lower()=='outfile10':
                                outtraj9 = value
                        elif option.lower()=='outfile11':
                                outtraj10 = value
                        elif option.lower()=='outfile12':
                                outtraj11 = value
                        elif option.lower()=='outfile13':
                                outtraj12 = value
                        elif option.lower()=='outfile14':
                                outtraj13 = value
                        else :
				print "Option:", option, " is not recognized"
        f.close()

# read in command line argument                                                                                                    
cfg_file = sys.argv[1]

# read cfg file                                                                                                                    
ParseConfigFile(cfg_file)
possibledimerresids = open(out1_file,'w')

# initiate MDAnalysis coordinate universe                                                                                          
u = MDAnalysis.Universe(top_file, traj_file)

#make selections between residues in dipeptide
NT = []
CT = []
CO = []
NH = []
PH = []
PH1 = []
PH2 = []
FF = []
for i in range(1,121, 2):
    NT.append("resid %s and name N H1 H2 H3" %(i))
    CO.append("resid %s and name C O" %(i))
    PH1.append("resid %s and name CG CD1 CD2 CE1 CE2 CZ" %(i))
for i in range(2,121, 2):
    CT.append("resid %s and name C OXT O" %(i))
    NH.append("resid %s and name N H" %(i))
    PH2.append("resid %s and name CG CD1 CD2 CE1 CE2 CZ" %(i))
for i in range(1,121, 2):
    j=i+1
    FF.append("resid %s or resid %s" %(i,j))
    #PH.append("(resid %s and name CB CG CD1 CD2 CE1 CE2 CZ) or (resid %s and name CB CG CD1 CD2 CE1 CE2 CZ)" %(i,j))

#cutoffs should be determined from histogramed pair distances from simulation
FFFFcutoff = 16.**2
NTNTcutoff = 3.8**2
NTCOcutoff = 4.6**2
NTNHcutoff = 5.0**2
NTCTcutoff = 4.3**2
COCOcutoff = 3.8**2
CONHcutoff = 5.8**2
COCTcutoff = 3.8**2
NHNHcutoff = 3.8**2
NHCTcutoff = 4.2**2
CTCTcutoff = 4.8**2
PHPHcutoff1 = 5.5**2
PHPHcutoff2 = 5.8**2
gencutoff = 4.0**2


FFdimer = u.select_atoms(FF[1],FF[2]) #needed so as to tell MDAnalysis the number of atoms for each "timestep"

adjacencymatrix = np.zeros((len(u.trajectory),1770)) 
count=-1
n=FFdimer.n_atoms
# Loop through trajectory 
with MDAnalysis.Writer(outtraj1, n) as Q, MDAnalysis.Writer(outtraj2, n) as R, MDAnalysis.Writer(outtraj3, n) as S, MDAnalysis.Writer(outtraj4, n) as T, MDAnalysis.Writer(outtraj5, n) as U, MDAnalysis.Writer(outtraj6, n) as V, MDAnalysis.Writer(outtraj7, n) as W, MDAnalysis.Writer(outtraj8, n) as X, MDAnalysis.Writer(outtraj9, n) as Y, MDAnalysis.Writer(outtraj10, n) as Z, MDAnalysis.Writer(outtraj11, n) as A, MDAnalysis.Writer(outtraj12, n) as B, MDAnalysis.Writer(outtraj13, n) as C:
    for ts in u.trajectory[::1]:
        count+=1
        temp=-1
        possibledimerresids.write('Frame %s of orginal trajectory\n' %(count))
        # calculate distance and find pairs of molecules within the cutoff distance
        for i in range(len(FF)):
            for j in range(len(FF)):
                if i < j:
                    temp+=1
                    FF1 = u.select_atoms(FF[i])
                    FF2 = u.select_atoms(FF[j])
                    comdist = computePbcDist2(FF1.center_of_mass(),FF2.center_of_mass(),u.dimensions[:3])   
                    if comdist <= FFFFcutoff:
                            NT1 = u.select_atoms(NT[i])
                            CO1 = u.select_atoms(CO[i])
                            NH1 = u.select_atoms(NH[i])
                            CT1 = u.select_atoms(CT[i])
                            NT2 = u.select_atoms(NT[j])
                            CO2 = u.select_atoms(CO[j])
                            NH2 = u.select_atoms(NH[j])
                            CT2 = u.select_atoms(CT[j])
                            dist2 = computePbcDist2(NT1.center_of_mass(),CO2.center_of_mass(),u.dimensions[:3])
                            dist4 = computePbcDist2(NT1.center_of_mass(),CT2.center_of_mass(),u.dimensions[:3])
                            dist5 = computePbcDist2(CO1.center_of_mass(),NT2.center_of_mass(),u.dimensions[:3])
                            dist7 = computePbcDist2(CO1.center_of_mass(),NH2.center_of_mass(),u.dimensions[:3])
                            dist10 = computePbcDist2(NH1.center_of_mass(),CO2.center_of_mass(),u.dimensions[:3])
                            dist12 = computePbcDist2(NH1.center_of_mass(),CT2.center_of_mass(),u.dimensions[:3])
                            dist13 = computePbcDist2(CT1.center_of_mass(),NT2.center_of_mass(),u.dimensions[:3])
                            dist15 = computePbcDist2(CT1.center_of_mass(),NH2.center_of_mass(),u.dimensions[:3])
                            if dist2 <= NTCOcutoff or dist4 <= NTCTcutoff or dist5 <= NTCOcutoff  or dist7 <= CONHcutoff or dist10<= CONHcutoff or dist12<= NHCTcutoff or dist13<= NTCTcutoff or dist15<= NHCTcutoff:
                                dist1 = computePbcDist2(NT1.center_of_mass(),NT2.center_of_mass(),u.dimensions[:3])
                                dist3 = computePbcDist2(NT1.center_of_mass(),NH2.center_of_mass(),u.dimensions[:3])
                                dist6 = computePbcDist2(CO1.center_of_mass(),CO2.center_of_mass(),u.dimensions[:3])
                                dist8 = computePbcDist2(CO1.center_of_mass(),CT2.center_of_mass(),u.dimensions[:3])
                                dist9 = computePbcDist2(NH1.center_of_mass(),NT2.center_of_mass(),u.dimensions[:3])
                                dist11 = computePbcDist2(NH1.center_of_mass(),NH2.center_of_mass(),u.dimensions[:3])
                                dist14 = computePbcDist2(CT1.center_of_mass(),CO2.center_of_mass(),u.dimensions[:3])
                                dist16 = computePbcDist2(CT1.center_of_mass(),CT2.center_of_mass(),u.dimensions[:3])
                                FFdimer1 = u.select_atoms(FF[i],FF[j])
                                Q.write(FFdimer1)
                                if dist4 <= NTCTcutoff and dist13 <= NTCTcutoff:
                                    possibledimerresids.write(' Anti-Aligned Stacked resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 2
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    R.write(FFdimer1)
                                elif dist1 <= NTNTcutoff and dist16 <= CTCTcutoff:
                                    possibledimerresids.write(' Aligned Stacked resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 3
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    A.write(FFdimer1)
                                elif (dist4 <= NTCTcutoff and dist3 >= NTNHcutoff and dist1 >= NTNTcutoff and dist2 >= NTCOcutoff and dist12 >= NHCTcutoff and dist15 >= NHCTcutoff) or (dist13 <= NTCTcutoff and dist9 >= NTNHcutoff and dist1 >= NTNTcutoff and dist2 >= NTCOcutoff and dist12 >= NHCTcutoff and dist15 >= NHCTcutoff):
                                    possibledimerresids.write(' Head-Tail resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 4
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    S.write(FFdimer1)
                                elif (dist2 <= NTCOcutoff and dist4 <= NTCTcutoff and dist15 >= NHCTcutoff and dist12 >= NHCTcutoff) or (dist5 <= NTCOcutoff and dist13 <= NTCTcutoff and dist12 >= NHCTcutoff and dist15 >= NHCTcutoff):
                                    possibledimerresids.write(' NT-AM T-shaped resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 5
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    T.write(FFdimer1)
                                elif (dist15 <= NHCTcutoff and dist13 <= NTCTcutoff and dist1 >= NTNTcutoff and dist5 >= NTCOcutoff and dist9 <= NTNHcutoff) or (dist12 <= NHCTcutoff and dist3 >= NTNHcutoff and dist4 <= NTCTcutoff and dist1 >= NTNTcutoff and dist2 >= NTCOcutoff):
                                    possibledimerresids.write(' CT-AM T-shaped resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 6
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    U.write(FFdimer1)
                                elif (dist4 <= NTCTcutoff and dist3 <= NTNHcutoff and dist15 >= NHCTcutoff and dist2 >= NTCOcutoff) or (dist13 <= NTCTcutoff and dist9 <= NTNHcutoff and dist12 >=NHCTcutoff and dist5 >= NTCOcutoff):
                                    possibledimerresids.write(' L-shaped resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 9
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    X.write(FFdimer1)
                                elif (dist15 <= NHCTcutoff and dist13 <= NTCTcutoff and dist5 <= NTCOcutoff) or (dist12 <= NHCTcutoff and dist4 <= NTCTcutoff and dist2 <= NTCOcutoff):
                                    possibledimerresids.write(' "3 Point" Parallel Displaced resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 7
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    V.write(FFdimer1)
                                elif (dist4 <= NTCTcutoff and dist12 <= NHCTcutoff) or (dist13 <= NTCTcutoff and dist15 <=NHCTcutoff):
                                #elif (dist4 <= NTCTcutoff and dist3 <= NTNHcutoff and dist12 <= NHCTcutoff) or (dist13 <= NTCTcutoff and dist9 <= NTNHcutoff and dist15 <=NHCTcutoff):
                                    possibledimerresids.write(' NT-CT "2 points" Parallel Displaced resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 8
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    W.write(FFdimer1)
                                elif dist2 <=NTCOcutoff or dist5 <=NTCOcutoff:
                                    possibledimerresids.write(' NT-CO resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 10
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    B.write(FFdimer1)
                                elif dist7 <=CONHcutoff or dist10 <=CONHcutoff:
                                    possibledimerresids.write(' NT-CO resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 11
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    C.write(FFdimer1)
                                else:
                                #elif (dist1 <= NTCTcutoff or dist2 <= NTCTcutoff or dist3 <= NTCOcutoff or dist4 <= NTCOcutoff or dist5 <= NHCTcutoff or dist6 <= NHCTcutoff):
                                    possibledimerresids.write(' REMAINDER NTNT %s NTCO %s NTNH %s NTCT %s CONT %s COCO %s CONH %s COCT %s NHNT %s NHCO %s NHNH %s NHCT %s CTNT %s CTCO %s CTNH %s CTCT %s\n ' %(np.sqrt(dist1), np.sqrt(dist2), np.sqrt(dist3), np.sqrt(dist4), np.sqrt(dist5), np.sqrt(dist6), np.sqrt(dist7), np.sqrt(dist8), np.sqrt(dist9), np.sqrt(dist10), np.sqrt(dist11), np.sqrt(dist12), np.sqrt(dist13), np.sqrt(dist14), np.sqrt(dist15), np.sqrt(dist16)))
                                    #possibledimerresids.write('REMAINDER resid %s %s %s %s\n ' %(i*2+1, i*2+2, j*2+1, (j+1)*2))
                                    adjacencymatrix[(count,temp)] = 12
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    Y.write(FFdimer1)
                            else:
                                PH11 = u.select_atoms(PH1[i])
                                PH21 = u.select_atoms(PH2[i])
                                PH12 = u.select_atoms(PH1[j])
                                PH22 = u.select_atoms(PH2[j])
                                dist17 = computePbcDist2(PH11.center_of_mass(),PH12.center_of_mass(),u.dimensions[:3])
                                dist18 = computePbcDist2(PH21.center_of_mass(),PH22.center_of_mass(),u.dimensions[:3])
                                dist19 = computePbcDist2(PH21.center_of_mass(),PH12.center_of_mass(),u.dimensions[:3])
                                dist20 = computePbcDist2(PH11.center_of_mass(),PH22.center_of_mass(),u.dimensions[:3])
                                if dist17 < PHPHcutoff1 or dist18 < PHPHcutoff2 or dist19 < PHPHcutoff2 or dist20 < PHPHcutoff2:
                                    adjacencymatrix[(count,temp)] = 1
                                    FFdimer1 = u.select_atoms(FF[i],FF[j])
                                    Z.write(FFdimer1)


Q.close()
R.close()
S.close()
T.close()
U.close()
V.close()
W.close()
X.close()
Y.close()
Z.close()
possibledimerresids.close()
np.savetxt("rep3adjacencymatrix.dat", adjacencymatrix, fmt='%d')
                        
