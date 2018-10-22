
import numpy as np
import MDAnalysis
import matplotlib.pyplot as plt
import cython

def gauss(x,sigma,xi):
    return np.exp(-0.5*(x-xi)**2/sigma**2)

def compute_pbc_r(coord1, coord2, box, hbox):
    r = np.empty(3,dtype=float) 
    for k in range(3):
        r[k] = coord1[k]-coord2[k]
        # apply periodic boundary conditions
        if r[k] > hbox[k]:
            r[k] -= box[k]
        elif r[k] < -hbox[k]:
            r[k] += box[k]
    return r

# main

E0AmideI = 1733.11*0.93                 # unperturbed energy of mode 
u0AmideI = np.array([0.81675418,0.57601332,0.03348519])
u0AmideI /= np.linalg.norm(u0AmideI)
dipoleStrengthAmideI = 0.374      # oscillator strength in Debye 

E0Cterm = 1664.14*0.93
u0Cterm = np.array([-0.6311285, 0.58748156, 0.50649998])
u0Cterm /= np.linalg.norm(u0Cterm)
dipoleStrengthCterm = 1.55*dipoleStrengthAmideI

E0AmideII = 1533.10*0.93                 # unperturbed energy of mode 
u0AmideII = np.array([-0.17023415, 0.97321288,-0.1545219 ])
u0AmideII /= np.linalg.norm(u0AmideII)
dipoleStrengthAmideII = 0.47*dipoleStrengthAmideI      # oscillator strength in Debye 

gamma = 10.0
hom = 5.0


# load trajectory using MDAnalysis
traj_file = "First250nsEvery25stepsOverTriplicate.dcd"
top_file = "60FF_ff14ipq_rep1.stripwater.prmtop"
coord = MDAnalysis.Universe(top_file, traj_file)
# make an atom selection
sel = []
sel.append("name C O")  # Amide C 0
sel.append("name N")  # Amide N
selCterm = "name C O OXT"

# define vector of dipoles
nAmide = coord.atoms.n_residues//2 # number of molecules of FF
nCterm = nAmide
N = 2*nAmide+nCterm
uAmideI = np.empty((nAmide,3),dtype=float)
dipoleCoordAmideI = np.empty((nAmide,3),dtype=float)
uAmideII = np.empty((nAmide,3),dtype=float)
dipoleCoordAmideII = np.empty((nAmide,3),dtype=float)
uCterm = np.empty((nCterm,3),dtype=float)
dipoleCoordCterm = np.empty((nCterm,3),dtype=float)
oscStrength = np.empty(N,dtype=float)
B = np.empty((N,N),dtype=float)
minX = 1400.0
maxX = 1800.0
deltaX = 0.01
x = np.arange(minX,maxX,deltaX)
f = np.zeros(x.size,dtype=float)
for ts in coord.trajectory:
    #if coord.trajectory.n_frames - ts.frame < 100:
        print(ts.frame, "of", coord.trajectory.n_frames)
        box = ts.dimensions[0:3]
        hbox = box/2.0
        # determine axes and dipole moment for each residue
        for resid1 in range(nAmide):
            # amide C and O are on odd residue
            selection = "(resid " + str(resid1*2+1) + " and " + sel[0] + ")"
            # amide N is on even residue
            selection += " or (resid " + str(resid1*2+2) + " and " + sel[1] + ")"
            currentAmide = coord.select_atoms(selection)
            # compute C=0 bond vector
            co = compute_pbc_r(currentAmide.positions[0,:],currentAmide.positions[1,:], box, hbox)
            # normalize
            coMag = np.linalg.norm(co)
            co /= coMag
            # compute C=N bond vector
            cn = compute_pbc_r(currentAmide.positions[0,:],currentAmide.positions[2,:], box, hbox)
            # normalize
            cn /= np.linalg.norm(cn)
            # x-axis is along co vector
            xp = np.copy(co)
            # z-axis is cross between C=O and C-N bond vectors
            zp = np.cross(co,cn)
            zp /= np.linalg.norm(zp)
            yp = np.cross(xp,zp)
            yp /= np.linalg.norm(yp)
            # place dipole in orientation of current molecule
            uAmideI[resid1,:] = u0AmideI[0]*xp + u0AmideI[1]*yp + u0AmideI[2]*zp
            uAmideII[resid1,:] = u0AmideII[0]*xp + u0AmideII[1]*yp + u0AmideII[2]*zp
            # save coordinate of Amide mode for current molecule
            dipoleCoordAmideI[resid1,:] = (currentAmide.positions[0,:] + currentAmide.positions[1,:])/2.0
            dipoleCoordAmideII[resid1,:] = currentAmide.positions[2,:]   # nitrogen position
            # place Cterm dipole in orientation of current molecule
            uCterm[resid1,:] = u0Cterm[0]*xp + u0Cterm[1]*yp + u0Cterm[2]*zp
            # save coordinate of cterm mode for current molecule
            selection = "resid " + str(resid1*2+2) + " and " + selCterm
            currentCterm = coord.select_atoms(selection)
            dipoleCoordCterm[resid1,:] = currentCterm.center_of_mass()

        uAmideI *= dipoleStrengthAmideI
        uAmideII *= dipoleStrengthAmideII
        uCterm *= dipoleStrengthCterm

        u = np.concatenate((uAmideI,uCterm,uAmideII),axis=0)
        dipoleCoord = np.concatenate((dipoleCoordAmideI,dipoleCoordCterm,dipoleCoordAmideII),axis=0)
        # calculate coupling
        rVec = np.empty(3,dtype=float)
        B[N-1,N-1] = 0.0  
        for i in range(N-1):
            B[i,i] = 0.0
            for j in range(i+1,N):
                # ignore intramolecular coupling
                if (j == i+nAmide or j == i+2*nAmide):
                    B[i,j] = 0.0
                else:    
                    rVec = dipoleCoord[i,:] - dipoleCoord[j,:]  # vector from center of TS on molecule j to molecule i
                    rMag = np.linalg.norm(rVec)
                    B[i,j] = np.dot(u[i,:],u[j,:])/(rMag**3) - 3.0*(np.dot(u[i,:],rVec)*np.dot(u[j,:],rVec))/(rMag**5) # transition dipole coupling
                    B[i,j] *= 5034.0  #? to convert to Debye 
                B[j,i] = B[i,j]   # symmetrize

        H0 = np.diag(np.concatenate((E0AmideI*np.ones(nAmide),E0Cterm*np.ones(nCterm),E0AmideII*np.ones(nAmide))))
        H1 = H0+B

        eigenValues, eigenVectors = np.linalg.eig(H1)

        u1 = np.matmul(u.T,eigenVectors).T
        oscStrength = np.sqrt(np.diag(np.matmul(u1,u1.T)))

        for i in range(N):
            f += oscStrength[i]*gauss(x,hom,eigenValues[i])


f /= np.amax(f)
np.savetxt("First250ns60ff_amide1_cterm_amide2_linear_params_cterm_com_irspectrum.txt",np.column_stack((x,f)))
