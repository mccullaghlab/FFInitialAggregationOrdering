import numpy as np

adjacencymatrix=np.loadtxt("3adjacencymatrix.dat")
clustersperframe=np.zeros((len(adjacencymatrix),60))

Nmono=60
for ts in range(len(adjacencymatrix)):
    nclu = 0
    clu=np.zeros(60,dtype=int)
    todo=np.zeros(60,dtype=int)
    for i in range(60):
        if clu[i] == 0:
            nclu+=1
            clu[i] = nclu
            todo[0] = i
            n=0
            m=0
            while n >= m:
                j = todo[m]
                for k in range(60):
                    if j <= k:
                        jp = k - 1 + ((j)*(2*Nmono-j-3))/2
                    else:
                        jp = j - 1 + ((k)*(2*Nmono-k-3))/2
                    #print ts, jp, j, k
                    #if adjacencymatrix[ts,jp] >= 3.0:
                    if adjacencymatrix[ts,jp] == 2.0 or adjacencymatrix[ts,jp] == 4.0:
                    #if adjacencymatrix[ts,jp] >= 2.0 and not adjacencymatrix[ts,jp] == 4.0:
                        if clu[k] == 0:
                            clu[k] = nclu
                            n+=1
                            todo[n]=k
                m+=1
    clustersperframe[ts]=clu
#np.savetxt("6clustersperframeevery25stepsGT3.dat",clustersperframe,fmt='%d')
np.savetxt("3clustersperframeevery25steps2and4.dat",clustersperframe,fmt='%d')
#np.savetxt("6clustersperframeevery25stepswithsidechains.dat",clustersperframe,fmt='%d')
