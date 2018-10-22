import os
import numpy as np
import MDAnalysis
import matplotlib.pyplot as plt
import matplotlib.markers as mmark


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
#####box may be an issue below
def computeOrthogVec(r1,r2,r3,box): 
        vec1=np.zeros(3, dtype=np.float)
        vec2=np.zeros(3, dtype=np.float)
        for i in range(0,3):
            vec1[i] = r2[i] - r1[i]
            vec2[i] = r3[i] - r2[i]
        crossvec = np.cross(vec1,vec2)
        crossvec /= np.linalg.norm(crossvec)
        return crossvec; 

# load trajectory using MDAnalysis
#traj_file = "Last250nsEvery25stepsRep2.dcd"
#top_file = "60FF_ff14ipq_rep1.stripwater.prmtop"
#coord = MDAnalysis.Universe(top_file, traj_file)
#
#
## make an atom selection
#PHEsel=["name CG","name CD1","name CD2","name CE1","name CE2","name CZ"]
#
#
#coord.atoms.n_residues
#
#
#
#
#firstCut = 16
#firstCut2 = firstCut*firstCut
#secondCut = 6
#secondCut2 = secondCut*secondCut
#count = 0
#dimerDataList = []
#dimerMetaDataList = []
#PHEang1 = []
#PHEang2 = []
#PHEang3 = []
#PHEang4 = []
#
#for ts in coord.trajectory[0::1]:
#        if ts.frame%1 == 0:
#            print(ts.frame, "of", coord.trajectory.n_frames)
#            box = ts.dimensions[0:3]
#            hbox = box/2.0
#            for resid1 in range(coord.atoms.n_residues/2-1):
#                selection1 = "resid " + str(resid1*2+1) + " " + str(resid1*2+2)
#                PHE1orthosel1 = "resid " + str(resid1*2+1) + " and name CE1"
#                PHE1orthosel2 = "resid " + str(resid1*2+1) + " and name CZ"
#                PHE1orthosel3 = "resid " + str(resid1*2+1) + " and name CE2"
#                PHE2orthosel1 = "resid " + str(resid1*2+2) + " and name CE1"
#                PHE2orthosel2 = "resid " + str(resid1*2+2) + " and name CZ"
#                PHE2orthosel3 = "resid " + str(resid1*2+2) + " and name CE2"
#                mol1 = coord.select_atoms(selection1)
#                for resid2 in range(resid1+1,coord.atoms.n_residues/2):
#                    selection2 = "resid " + str(resid2*2+1)+ " " + str(resid2*2+2)
#                    PHE3orthosel1 = "resid " + str(resid2*2+1) + " and name CE1"
#                    PHE3orthosel2 = "resid " + str(resid2*2+1) + " and name CZ"
#                    PHE3orthosel3 = "resid " + str(resid2*2+1) + " and name CE2"
#                    PHE4orthosel1 = "resid " + str(resid2*2+2) + " and name CE1"
#                    PHE4orthosel2 = "resid " + str(resid2*2+2) + " and name CZ"
#                    PHE4orthosel3 = "resid " + str(resid2*2+2) + " and name CE2"
#                    mol2 = coord.select_atoms(selection2)
#                    dist2 = compute_pbc_dist2(mol1.center_of_geometry(),mol2.center_of_geometry(),box,hbox)
#                    if dist2 < firstCut2:  # we might have a dimer
#                        tempList = []
#                        orthotempList = []
#                        for i in range(len(PHEsel)):
#                            if i < len(PHEsel)/2:
#                                selection = "resid " + str(resid1*2+1) + " and " + PHEsel[i]
#                            else:
#                                selection = "resid " + str(resid1*2+2) + " and " + PHEsel[i]
#                            subMol1 = coord.select_atoms(selection)
#                            for j in range(len(PHEsel)):
#                                if j < len(PHEsel)/2:
#                                    selection = "resid " + str(resid2*2+1) + " and " + PHEsel[j]
#                                else:
#                                    selection = "resid " + str(resid2*2+2) + " and " + PHEsel[j]
#                                subMol2 = coord.select_atoms(selection)
#                                dist2 = compute_pbc_dist2(subMol1.center_of_mass(), subMol2.center_of_mass(), box, hbox)
#                                tempList.append(dist2)
#                        if min(tempList) < secondCut2:  # now we are sure we have a dimer
#                            #compute orthogonal vector of each phenyl ring in the dimer (4 vectors)
#                            vec1 = computeOrthogVec(coord.select_atoms(PHE1orthosel1).center_of_mass(),coord.select_atoms(PHE1orthosel2).center_of_mass(),coord.select_atoms(PHE1orthosel3).center_of_mass(),box)
#                            vec2 = computeOrthogVec(coord.select_atoms(PHE2orthosel1).center_of_mass(),coord.select_atoms(PHE2orthosel2).center_of_mass(),coord.select_atoms(PHE2orthosel3).center_of_mass(),box)
#                            vec3 = computeOrthogVec(coord.select_atoms(PHE3orthosel1).center_of_mass(),coord.select_atoms(PHE3orthosel2).center_of_mass(),coord.select_atoms(PHE3orthosel3).center_of_mass(),box)
#                            vec4 = computeOrthogVec(coord.select_atoms(PHE4orthosel1).center_of_mass(),coord.select_atoms(PHE4orthosel2).center_of_mass(),coord.select_atoms(PHE4orthosel3).center_of_mass(),box)
#                            #compute angle between vectors of phenyl rings
#                            vecsdot1 = np.dot(vec1,vec3)
#                            vecsdot2 = np.dot(vec1,vec4)
#                            vecsdot3 = np.dot(vec2,vec3)
#                            vecsdot4 = np.dot(vec2,vec4)
#                            
#                            if vecsdot1 >= 1.0:
#                                PHEang1.append(0.0)
#                            elif vecsdot1 <= -1.0:
#                                PHEang1.append(180.)
#                            else:
#                                PHEang1.append(np.clip(np.arccos(vecsdot1)*180./np.pi,0.,180.))
#                                #ang= np.arccos(vecsdot1)
#                                #ang2 = ang/np.sin(ang)
#                                #PHEang1.append(ang2*180./np.pi)
#                            if vecsdot2 >= 1.0:
#                                PHEang2.append(0.0)
#                            elif vecsdot2 <= -1.0:
#                                PHEang2.append(180.)
#                            else:
#                                PHEang2.append(np.arccos(vecsdot2)*180./np.pi)
#                            if vecsdot3 >= 1.0:
#                                PHEang3.append(0.0)
#                            elif vecsdot3 <= -1.0:
#                                PHEang3.append(180.)
#                            else:
#                                PHEang3.append(np.arccos(vecsdot3)*180./np.pi)
#                            if vecsdot4 >= 1.0:
#                                PHEang4.append(0.0)
#                            elif vecsdot4 <= -1.0:
#                                PHEang4.append(180.)
#                            else:
#                                #ang= np.arccos(vecsdot4)
#                                #ang2 = ang/np.sin(ang)
#                                #PHEang4.append(ang2*180./np.pi)
#                                PHEang4.append(np.arccos(vecsdot4)*180./np.pi)
#
#                            #print vecsdot4, PHEang4, selection1, selection2
#                            
#                            # store dimer distance data
#                            #dimerDataList.append([])
#                            #for el in tempList:
#                            #    dimerDataList[count].append(np.sqrt(el))
#                            # store dimer meta data
#                            dimerMetaDataList.append([])
#                            dimerMetaDataList[count].append(ts.frame)
#                            dimerMetaDataList[count].append(selection1)
#                            dimerMetaDataList[count].append(selection2)
#                            count += 1
#                
#                
#for i in range(len(PHEang4)):
#    if PHEang4[i] >90.0:
#        PHEang4[i]=90.0 - 90.*(PHEang4[i]/90.-1.)
#for i in range(len(PHEang1)):
#    if PHEang1[i] >90.0:
#        PHEang1[i]=90.0 - 90.*(PHEang1[i]/90.-1.)
#
#
#
#PHEData1 = np.asmatrix(PHEang1)
#PHEData2 = np.asmatrix(PHEang2)
#PHEData3 = np.asmatrix(PHEang3)
#PHEData4 = np.asmatrix(PHEang4)
#np.savetxt("PHEang1.60FFLast250ns_rep2.dat",PHEData1)
#np.savetxt("PHEang2.60FFLast250ns_rep2.dat",PHEData2)
#np.savetxt("PHEang3.60FFLast250ns_rep2.dat",PHEData3)
#np.savetxt("PHEang4.60FFLast250ns_rep2.dat",PHEData4)
#metaDataOut = open("PHEangs.60FFLast250ns_rep2.metadata","w")
#for i in range(count):
#    metaDataOut.write("%d %s %s\n" % (dimerMetaDataList[i][0], dimerMetaDataList[i][1],dimerMetaDataList[i][2]))
#metaDataOut.close()
## read data if already produced
#PHEData1 = np.loadtxt("PHEang1.FFLast250.dat")
#PHEData2 = np.loadtxt("PHEang2.FFLast250.dat")
#PHEData3 = np.loadtxt("PHEang3.FFLast250.dat")
#PHEData4 = np.loadtxt("PHEang4.FFLast250.dat")
PHEData1 = np.loadtxt("PHEang1.60FFLast250ns_rep3.dat")
PHEData2 = np.loadtxt("PHEang2.60FFLast250ns_rep3.dat")
PHEData3 = np.loadtxt("PHEang3.60FFLast250ns_rep3.dat")
PHEData4 = np.loadtxt("PHEang4.60FFLast250ns_rep3.dat")
for i in range(len(PHEData2)):
    if PHEData2[i] >90.0:
        PHEData2[i]=180.0 - PHEData2[i]

for i in range(len(PHEData1)):
    PHEData1[i]= float(PHEData1[i]*np.sin(np.deg2rad(PHEData1[i])))
for i in range(len(PHEData2)):
    PHEData2[i]= float(PHEData2[i]*np.sin(np.deg2rad(PHEData2[i])))
for i in range(len(PHEData4)):
    PHEData4[i]= float(PHEData4[i]*np.sin(np.deg2rad(PHEData4[i])))

plt.figure(figsize=(15,15))
plt.rcParams['axes.linewidth'] = 3.7
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
plt.rc('axes', axisbelow=True)
handles1 = plt.Rectangle((0,0),1,1,fc="b")
handles2 = plt.Rectangle((0,0),1,1,fc="g")
handles3 = plt.Rectangle((0,0),1,1,fc="k")
labels= ["NTPHE-NTPHE","NTPHE-CTPHE","CTPHE-CTPHE"]
plt.hist(PHEData1, normed=1, bins=180, color='b',histtype='step', linewidth=4)
plt.hist(PHEData2, normed=1, bins=180, color='g',histtype='step', linewidth=4)
#plt.hist(PHEData3.A1, normed='True', bins=360, color='r',histtype='step', linewidth=3)
plt.hist(PHEData4, normed=1, bins=180, color='k',histtype='step', linewidth=4)
plt.legend([handles1, handles2, handles3],labels, fontsize='31', loc='upper right', ncol=1, numpoints = 0.1)
#plt.legend(['NTPHE-NTPHE','CTPHE-CTPHE'], fontsize='32', loc='upper right', ncol=1, numpoints = 0.5)
#plt.legend(['NTPHE-NTPHE','NTPHE-CTPHE','CTPHE-NTPHE','CTPHE-CTPHE'], fontsize='20', loc='upper right', ncol=1, numpoints = 0.5)
#plt.title(r'%s' %(angle), size='30')
plt.xlabel('Dihedral ($\degree$)', size='40')
plt.ylabel('Probability density',size='40')
plt.tick_params(labelsize=40)
plt.xlim([0,90])
plt.ylim([0,.05])
plt.xticks(range(0,100,30))
#plt.yticks(np.arange(0,.051,step=0.004))
plt.savefig('PHEangs_60FFLast250rep3.png',transparent='True')#,dpi=600)
plt.savefig('PHEangs_60FFLast250rep3.eps',format='eps', bbox_inches='tight')
plt.savefig('PHEangs_60FFLast250rep3.pdf')
plt.close()     


