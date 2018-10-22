import numpy as np
import MDAnalysis
import time
import math
import sys

# main program
def main():
    start = time.time()
    coord1 = MDAnalysis.Universe("2FF_ff14ipq.prmtop","w07/2FF_ff14ipq.w07.run00.C.nc",forces=True)
    coord2 = MDAnalysis.Universe("2FF_ff14ipq.prmtop","w07/2FF_ff14ipq.w07.run00.LJ.nc",forces=True)
    coord3 = MDAnalysis.Universe("2FF_ff14ipq.stripwater.prmtop","w07/2FF_ff14ipq.w07.run00_nowater.C.nc",forces=True)
    coord4 = MDAnalysis.Universe("2FF_ff14ipq.stripwater.prmtop","w07/2FF_ff14ipq.w07.run00_nowater.LJ.nc",forces=True)
    solute1 = [coord1.select_atoms("resid 1 2"), coord1.select_atoms("resid 3 4")]
    solute2 = [coord2.select_atoms("resid 1 2"), coord2.select_atoms("resid 3 4")]
    solute3 = [coord3.select_atoms("resid 1 2"), coord3.select_atoms("resid 3 4")]
    solute4 = [coord4.select_atoms("resid 1 2"), coord4.select_atoms("resid 3 4")]
    r=np.zeros((2,3))
    fCtot=np.zeros((2,3))
    fLJtot=np.zeros((2,3))
    fC=np.zeros((2,3))
    fLJ=np.zeros((2,3))
    #
    binSize=0.1
    rmax=17.5
    rmax2=rmax*rmax
    numBin=int(rmax/binSize)
    nHist=np.zeros(numBin)
    fCsolv=np.zeros(numBin)
    fLJsolv=np.zeros(numBin)
    fCsolu=np.zeros(numBin)
    fLJsolu=np.zeros(numBin)
    fCsolv2=np.zeros(numBin)
    fLJsolv2=np.zeros(numBin)
    fCsolu2=np.zeros(numBin)
    fLJsolu2=np.zeros(numBin)
    dims=np.zeros(3)
    hdims=np.zeros(3)
    #
    rlist=["03.50","04.00","04.50","05.00","05.50","06.00","06.50","07.00","07.50","08.00","08.50","09.00","09.50","10.00","10.50","11.00","11.50","12.00","12.50","13.00","13.50","14.00","14.50","15.00","15.50","16.50","17.00","17.50"]
    ilist=["00","01","02","03","04"]
    wlist=["07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35"]
    nr=0
    for w0 in wlist:
        for ip in ilist:
            nr+=1
            coord1.load_new("w"+w0+"/2FF_ff14ipq.w"+w0+".run"+ip+".C.nc",forces=True)
            coord2.load_new("w"+w0+"/2FF_ff14ipq.w"+w0+".run"+ip+".LJ.nc",forces=True)
            coord3.load_new("w"+w0+"/2FF_ff14ipq.w"+w0+".run"+ip+"_nowater.C.nc",forces=True)
            coord4.load_new("w"+w0+"/2FF_ff14ipq.w"+w0+".run"+ip+"_nowater.LJ.nc",forces=True)
            for ts in coord1.trajectory:
                ts2=coord2.trajectory[ts.frame-1]
                ts3=coord3.trajectory[ts.frame-1]
                ts4=coord4.trajectory[ts.frame-1]
                dims=coord1.dimensions
                hdims=dims/2.
                sys.stdout.write("Progress: {1}/{2}  {0:.2f}% Complete\r".format((float(ts.frame) / float(len(coord1.trajectory))) * 100,nr,len(wlist)*len(ilist)))
                sys.stdout.flush()
                for i in range(2):
                    r[i]=np.zeros(3)
                    mtot=0.
                    for b in solute1[i].atoms:
                        r[i]+=b.mass*b.position
                        mtot+=b.mass
                    r[i]/=mtot
                    fCtot[i]=np.zeros(3)
                    for b in solute1[i].atoms:
                        fCtot[i]+=b.force
                    fLJtot[i]=np.zeros(3)
                    for b in solute2[i].atoms:
                        fLJtot[i]+=b.force
                    fC[i]=np.zeros(3)
                    for b in solute3[i].atoms:
                        fC[i]+=b.force
                    fLJ[i]=np.zeros(3)
                    for b in solute4[i].atoms:
                        fLJ[i]+=b.force
                r12=r[0]-r[1]
                for i in range(3):
                    while r12[i] > hdims[i]:
                        r12[i]-=dims[i]
                    while r12[i] < -hdims[i]:
                        r12[i]+=dims[i]
                dist=np.dot(r12,r12)
                if dist < rmax2:
                    dist=np.sqrt(dist)
                    r12/=dist
                    ibin=int(dist/binSize)
                    nHist[ibin]+=2
                    fnow=np.dot(fCtot[0]-fC[0],r12)
                    fCsolv[ibin]+=fnow
                    fCsolv2[ibin]+=fnow*fnow
                    fnow=np.dot(fLJtot[0]-fLJ[0],r12)
                    fLJsolv[ibin]+=fnow
                    fLJsolv2[ibin]+=fnow*fnow
                    fnow=np.dot(fC[0],r12)
                    fCsolu[ibin]+=fnow
                    fCsolu2[ibin]+=fnow*fnow
                    fnow=np.dot(fLJ[0],r12)
                    fLJsolu[ibin]+=fnow
                    fLJsolu2[ibin]+=fnow*fnow
                    fnow=np.dot(fCtot[1]-fC[1],r12)
                    fCsolv[ibin]-=fnow
                    fCsolv2[ibin]+=fnow*fnow
                    fnow=np.dot(fLJtot[1]-fLJ[1],r12)
                    fLJsolv[ibin]-=fnow
                    fLJsolv2[ibin]+=fnow*fnow
                    fnow=np.dot(fC[1],r12)
                    fCsolu[ibin]-=fnow
                    fCsolu2[ibin]+=fnow*fnow
                    fnow=np.dot(fLJ[1],r12)
                    fLJsolu[ibin]-=fnow
                    fLJsolu2[ibin]+=fnow*fnow
    outFile = open("2FF_ff14ipq_frc.dat", 'w')
    outFile.write("   x      Coulomb Solvent   Variance C Solvent     LJ Solvent         LJ Variance      Coulomb Solute   Variance C Solute      LJ Solute         LJ Variance           Hist\n")
    for i in range(numBin):
        x=(i+0.5)*binSize
        if nHist[i] > 0.5:
            fCsolv[i] /=nHist[i]
            fLJsolv[i]/=nHist[i]
            fCsolu[i] /=nHist[i]
            fLJsolu[i]/=nHist[i]
            fCsolv2[i] /=nHist[i]
            fLJsolv2[i]/=nHist[i]
            fCsolu2[i] /=nHist[i]
            fLJsolu2[i]/=nHist[i]
            fCsolv2[i] =np.sqrt((fCsolv2[i] -fCsolv[i] *fCsolv[i] )/nHist[i])
            fLJsolv2[i]=np.sqrt((fLJsolv2[i]-fLJsolv[i]*fLJsolv[i])/nHist[i])
            fCsolu2[i] =np.sqrt((fCsolu2[i] -fCsolu[i] *fCsolu[i] )/nHist[i])
            fLJsolu2[i]=np.sqrt((fLJsolu2[i]-fLJsolu[i]*fLJsolu[i])/nHist[i])
        outFile.write("{:6.3f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:18.12f} {:12d}\n".format(x,fCsolv[i],fCsolv2[i],fLJsolv[i],fLJsolv2[i],fCsolu[i],fCsolu2[i],fLJsolu[i],fLJsolu2[i],int(nHist[i])))
    outFile.close()

    end = time.time()
    t = end - start
    print "\nTotal running time: {:.2f} sec".format(t)

# Main program code
main()
