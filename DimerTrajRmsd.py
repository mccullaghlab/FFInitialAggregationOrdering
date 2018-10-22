import numpy as np
from scipy.linalg import logm, expm
import MDAnalysis
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rms import rmsd

###README###
##I loaded the dimer type trajectory into vmd and then wrapped the traj. around resid 1 2. I then used the rmsd traj tool to align the traj around resid 1 2 backbone. I then wrote the centered trajectory and used this new traj. to find the rmsd. 


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


top_file = "PossibleDimers.prmtop"
#traj_file = "3ParallelStackedDimerscentered.dcd"
#traj_file = "3NTAMTshapedDimerscentered.dcd"
#traj_file = "3CTAMTshapedDimerscentered.dcd"
#traj_file = "3StaggeredW3pointsofcontactcentered.dcd"
#traj_file = "3StaggeredW2pointsofcontactcentered.dcd"
#traj_file = "3NTCOtshapedcentered.dcd"
#traj_file = "3CTNHtshapedcentered.dcd"
#traj_file = "3LshapedDimerscentered.dcd"
#traj_file = "3HeadTailDimerscentered.dcd"
traj_file = "3Remaindercentered.dcd"
#ref_pdb = "Refparallelstackeddimer.pdb"
#ref_pdb = "RefNTAMTshaped.pdb"
#ref_pdb = "RefCTAMTshaped.pdb"
#ref_pdb = "RefStaggeredW3points.pdb"
#ref_pdb = "RefStaggeredW2points.pdb"
#ref_pdb = "RefNTCOtshaped.pdb"
#ref_pdb = "RefCTNHtshaped.pdb"
#ref_pdb = "RefLshaped.pdb"
#ref_pdb = "RefHeadTail.pdb"
ref_pdb = "RefRemainder.pdb"
ref = MDAnalysis.Universe(ref_pdb)
#u = MDAnalysis.Universe(traj_file)
u = MDAnalysis.Universe(top_file, traj_file)

dimerRmsdDataList= np.zeros(len(u.trajectory))
count=0
for ts in u.trajectory:
    monomer1 = u.select_atoms("resid 1 and backbone or (resid 2 and name N CA C)")
    monomer2 = u.select_atoms("resid 3 and backbone or (resid 4 and name N CA C)")
    dimer1_positions=np.concatenate((np.asarray(monomer1.positions),np.asarray(monomer2.positions)),axis=0)
    dimer2_positions=np.concatenate((np.asarray(monomer2.positions),np.asarray(monomer1.positions)),axis=0)
    ref1_mon= ref.select_atoms("resid 1 and backbone or (resid 2 and name N CA C)")
    ref2_mon= ref.select_atoms("resid 3 and backbone or (resid 4 and name N CA C)")
    ref_positions= np.concatenate((np.asarray(ref1_mon.positions),np.asarray(ref2_mon.positions)),axis=0)
    firstrmsd = rmsd(dimer1_positions, ref_positions,superposition=True)
    secondrmsd = rmsd(dimer2_positions, ref_positions,superposition=True)
    dimerRmsdDataList[count] = min(firstrmsd,secondrmsd)
    count+=1
print "AVG:", np.mean(dimerRmsdDataList)
np.savetxt("remainderRMSDs.dat",dimerRmsdDataList)
#np.savetxt("headtailRMSDs.dat",dimerRmsdDataList)
#np.savetxt("lshapedRMSDs.dat",dimerRmsdDataList)
#np.savetxt("ctnhtshapedRMSDs.dat",dimerRmsdDataList)
#np.savetxt("ntcotshapedRMSDs.dat",dimerRmsdDataList)
#np.savetxt("staggeredw2pointsRMSDs.dat",dimerRmsdDataList)
#np.savetxt("ctamTshapedRMSDs.dat",dimerRmsdDataList)
#np.savetxt("nt_ct_stackedRMSDs.dat",dimerRmsdDataList)
