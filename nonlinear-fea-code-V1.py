# -*- coding: utf-8 -*-
"""
Spyder Editor

Change in V2 is addition of Displacement controlled BC
"""
#sum = 0
#n = int(input('Enter the value '))
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.tri as tri
from mesh import ReadAbaqusInp
import time
startTime = time.time()





#---------------------functions-------------------------
def Bmatrix(elabel, ElemData, NodeData, centroids):
    elem = ElemData[elabel]
    xi, xj, xm = NodeData[elem[0]-1,1], NodeData[elem[1]-1,1], NodeData[elem[2]-1,1]#, 0, 20, 0
    yi, yj, ym = NodeData[elem[0]-1,2], NodeData[elem[1]-1,2], NodeData[elem[2]-1,2]
    Adet = np.array([[1, xi, yi], [1, xj, yj], [1, xm, ym]])
    A = (1/2)*np.linalg.det(Adet)
    beta_i =  yj-ym 
    beta_j = ym-yi 
    beta_m = yi-yj
    gamma_i = -(xj-xm)
    gamma_j = -(xm-xi)
    gamma_m = -(xi-xj)
    b1 = [beta_i, 0, beta_j, 0, beta_m, 0]
    b2 = [0, gamma_i, 0, gamma_j, 0, gamma_m]
    b3 = [gamma_i, beta_i, gamma_j, beta_j, gamma_m, beta_m]
    Bmat = (1/(2*A))*np.array([b1, b2, b3])
    centroids[elabel, 0], centroids[elabel, 1] = (xi+xj+xm)/3, (yi+yj+ym)/3
    #centroids[elabel, 0]
    return Bmat, A
def TangentMatrix(E, nu):
    
    #D = E/(1-nu**2)*np.array([[1, nu, 0],[nu, 1, 0],[0, 0, (1-nu)/2]]) for plane stress
    D = E/((1-2*nu)*(1+nu))*np.array([[1-nu, nu, 0],[nu, 1-nu, 0],[0, 0, (1-2*nu)/2]])# for plane strain
    return D
def GlobalAssembly(ElemData, NodeData, LoadNode, E, nu, rglobal, kglobal):
    """Assembly of tri element stiffness (ke) and nodal force vector (rlocal) into 
    the global form rglobl, kglobal respectively.
    rlocal: elemental nodal force vector
    ke: elemental stiffness matrix 
    kglobal: global stiffness matrix 
    """
    for el in range(nele):
        elem = ElemData[el]
        B, A = Bmatrix(el, ElemData, NodeData, centroids)
        D = TangentMatrix(E, nu)
        tempM = np.dot(np.transpose(B),D)
        ke = A*np.dot(tempM,B) 
        connect = np.repeat(elem, 2)
        r = np.array([2, 1, 2, 1, 2, 1])
        rlocal = list(LoadNode[elem[0]-1,:]) + list(LoadNode[elem[1]-1,:])+ list(LoadNode[elem[2]-1,:])
        #print(rlocal)
        for i in range(0, elemdof):
            rglobal[(2*connect[i]-r[i]), 0] = rglobal[(2*connect[i]-r[i]), 0] + rlocal[i]
            #rglobal[(2*connect[i]-r[i]), 0] = rlocal[i]
            #print(rglobal)
            #print('\n')
            for j in range(0,elemdof):
                kglobal[(2*connect[i]-r[i]), (2*connect[j]-r[j])] = kglobal[(2*connect[i]-r[i]), (2*connect[j]-r[j])] + \
                ke[i, j]
                #print(kglobal[(2*connect[i]-r[i]), (2*connect[j]-r[j])])
                #print(rglobal)
                #print('\n')
    return kglobal, rglobal 
def DeAssemblyGlobalMatrix(kglobal, rglobal, LoadBCnodes, SupportBCnodes):
    kreduced, rreduced = kglobal.copy(), rglobal.copy()
    print('Displacement boundary condition successfully applied!')
    rreduced = np.ma.asarray(rreduced)
    kreduced = np.ma.asarray(kreduced)
    
    for icount in LoadBCnodes:
        rreduced = rreduced - xdisp*kreduced[:, icount].reshape(len(kreduced[:, icount]), 1)
        rreduced[icount, 0] = np.ma.masked
        kreduced[icount, icount] = np.ma.masked
    for icount in SupportBCnodes:
        rreduced[icount, 0] = np.ma.masked
        kreduced[icount, icount] = np.ma.masked 
    kreduced = np.ma.mask_rowcols(kreduced)
    rreduced =rreduced.compressed()
    kreduced =kreduced.compressed()
    kreduced = kreduced.reshape(len(rreduced), len(rreduced))
    
    return kreduced, rreduced
def SolveU(kreduced, rreduced, UnknownBCnodes, LoadBCnodes, SupportBCnodes, DisplacewithBC, xdisp):
    unknownU = np.linalg.solve(kreduced,rreduced)
    Displacement = DisplacewithBC.copy()
    print('Solution completed for unknown displacements')
    for indx in range(len(UnknownBCnodes)):
        index = int(UnknownBCnodes[indx])
        Displacement[index] =unknownU[indx] 
    for indx in range(len(LoadBCnodes)):
        index = int(LoadBCnodes[indx])
        Displacement[index] = xdisp#unknownU[indx] 
    for indx in range(len(SupportBCnodes)):
        index = int(SupportBCnodes[indx])
        Displacement[index] = 0#unknownU[indx] 
    return Displacement
def ComputeStressStrain(ElemData, nele, NodeData, Dispx, Dispy, E, nu):
    strain, stress = np.zeros([nele, 4]), np.zeros([nele, 4])
    for e in range(len(ElemData)):
        enum = ElemData[e]-1
        # elemDisp = [enum[0], enum[1], enum[2]]   
        B, Area = Bmatrix(e, ElemData, NodeData, centroids)     
        uarray = np.transpose(np.array([Dispx[enum[0]], Dispy[enum[0]],\
                                        Dispx[enum[1]], Dispy[enum[1]],\
                                        Dispx[enum[2]], Dispy[enum[2]]]))
        strain[e, 0:3] = np.matmul(B, uarray)   
        D = TangentMatrix(E, nu)
        stress[e, 0:3] = np.matmul(D, strain[e, 0:3])
        strain[e, 3] = 0
        stress[e, 3] = E*nu/((1-2*nu)*(1+nu))*(strain[e, 0] + strain[e, 1])
    return stress, strain    




###########################################################
#-------------------------------------------------------
TypeBC = 'Displacement'
#-------------------------------------------------------
nu = 0.3
E = 300 * 10**3
#NodeData = np.array([[1, 0, 0], [2, 0, 10], \
             #3, 20, 10], [4, 20, 0]])
print('Input file successfully imported')         
print('Reading geometry and mesh data!')
NodeData, ElemDatawithElemNo = ReadAbaqusInp('Job-11.inp')
#NodeData, ElemDatawithElemNo = MeshComputationfromPy()
#ElemDatawithElemNo = np.loadtxt('Elem-221.txt', delimiter = ',')
ElemData  = ElemDatawithElemNo[:,1:4]
ElemData = ElemData.astype(int)
centroids = np.zeros([len(ElemData), 2])

NodeCoord = NodeData[:, 1:3]
NodeintLabel = NodeData[:,0].astype(int)-1
ElemintLabel = ElemDatawithElemNo[:,0]
#plt.triplot(NodeCoord[:,0], NodeCoord[:,1], ElemData-1)
#--------------------------------------------------------
nele = len(ElemData)
ndof = 2
nnode = len(NodeData)
elemdof =6
LoadNode = np.zeros([nnode, 2])
rglobal  = np.zeros([nnode*ndof, 1])
kglobal = np.zeros([nnode*ndof,nnode*ndof])
print('Total no. of nodes : ', nnode)
print('Total no. of elements : ', nele)
print('Applying force on the nodes')

if TypeBC=='Force':
    for i in range(len(NodeCoord)):
        if NodeCoord[i, 0] == 20.0:
            LoadNode[i, 0] =200


print('Elementwise solution started: Building global force and stiffness matrix')
################################################################################
DisplacewithBC = np.zeros([ndof*nnode,1])
xdisp = 30 # prescribed displacement BC

print('Applying displacement boundary condition to the constrained node')
whereBC = [20, 0]
def getNodeset(whereBC, DisplacewithBC):
    SupportBCnodes = np.array([]).astype(int)
    LoadBCnodes = np.array([]).astype(int)
    UnknownBCnodes = np.array([]).astype(int)
    for i in NodeintLabel:
        if NodeCoord[i, 0] ==0.0:#NodeCoord[i, 0] <= 5.0 and NodeCoord[i, 1] == 0.0 : #and 
            DisplacewithBC[2*i, 0] = 0
            DisplacewithBC[2*i+1, 0] = 0
            SupportBCnodes = np.append(SupportBCnodes, [2*i, 2*i+1])
   # if NodeCoord[i, 1] ==0.0:
    #    DisplacewithBC[2*i+1, 0] = 0
    #if TypeBC == 'Displacement':
        if NodeCoord[i, 0] == whereBC[0]: #and NodeCoord[i, 1] == 0.0 :
            DisplacewithBC[2*i, 0] = 1
            print(NodeCoord[i, 1])
            #DisplacewithBC[2*i+1, 0] = xdisp
            LoadBCnodes = np.append(LoadBCnodes, 2*i)
#else:
        UnknownBCnodes = np.append(UnknownBCnodes, [2*i, 2*i+1])
    UnknownBCnodes = np.setxor1d(UnknownBCnodes, SupportBCnodes)
    UnknownBCnodes = np.setxor1d(UnknownBCnodes, LoadBCnodes)
    return SupportBCnodes, LoadBCnodes, UnknownBCnodes, DisplacewithBC
SupportBCnodes, LoadBCnodes, UnknownBCnodes, DisplacewithBC =  getNodeset(whereBC, DisplacewithBC)  

    
def InitalDisp(DisplacewithBC):
    DispVec = np.linspace(0, 10, 10)    
    DispMat = DispVec*DisplacewithBC   
    return DispMat
  #  if NodeCoord[i, 1] ==10.0:
   #     DisplacewithBC[2*i+1, 0] = 0
DispMat = InitalDisp(DisplacewithBC)
def DispInc(DispMat, ninc):
    du = DispMat[:,ninc+1] - DispMat[:, ninc]
    return du








################################################################################



kglobal, rglobal = GlobalAssembly(ElemData, NodeData, LoadNode, E, nu, rglobal, kglobal)
    
print('Global stiffness matrix successfully created')  


kreduced, rreduced = DeAssemblyGlobalMatrix(kglobal, rglobal, LoadBCnodes, SupportBCnodes)

print('Solving for unknown nodal displacements')  

Displacement = SolveU(kreduced, rreduced, UnknownBCnodes, LoadBCnodes, SupportBCnodes, DisplacewithBC, xdisp)
force = np.matmul(kglobal, Displacement)
#Displacement = np.concatenate((knownU, unknownU), axis = None)
#print('Global Matrix : \n', Displacement)
#print(DisplacewithBC)
DeformedCoord = np.zeros([nnode, 2])
print('Obtaining the deformed shape')

for ii in range(nnode):
    DeformedCoord[ii, :] = NodeCoord[ii] + np.array([Displacement[2*ii, 0], Displacement[2*ii+1, 0]])
    
DeformNodeData= np.zeros([len(DeformedCoord), 3])
DeformNodeData[:, 0] = np.linspace(1, len(DeformedCoord), len(DeformedCoord))-1
DeformNodeData[:, 1] = DeformedCoord[:, 0]
DeformNodeData[:, 2] = DeformedCoord[:, 1]
Dispx, forceX = np.array([]), np.array([])
Dispy, forceY = np.array([]), np.array([])
for i in range(nnode):
    Dispx = np.append(Dispx, Displacement[2*i])
    Dispy = np.append(Dispy, Displacement[2*i+1])
    forceX = np.append(forceX, force[2*i])
    forceY = np.append(forceY, force[2*i+1])
    
    
print('Calculating stress-strain response from displacement!')



stress, strain =  ComputeStressStrain(ElemData, nele, NodeData, Dispx, Dispy, E, nu)    
# Plotting the deformed geometry  
#MeshPlot(NodeCoord, ElemData, 'red', '-')
print('Analysis successfully completed!')
#MeshPlot(NodeCoord, ElemData, 'red', '-')
#MeshPlot(DeformedCoord, ElemData, 'blue', '--')
#fig, ax = plt.subplots()
#xcoord = NodeCoord[:,0]
#ycoord = NodeCoord[:,1]
deformcentroids = np.zeros([len(ElemData), 2])
for e in range(len(ElemData)):
        enum = ElemData[e]
        # elemDisp = [enum[0], enum[1], enum[2]]   
        B, Area = Bmatrix(e, ElemData, DeformNodeData, deformcentroids)  



contr1 = plt.triplot(NodeCoord[:,0], NodeCoord[:,1], ElemData-1,\
           'k-',lw = 0.30)
#plt.scatter(centroids[:, 0], centroids[:,1])
#x, y = np.meshgrid(centroids[:, 0], centroids[:, 1])
#plt.scatter(x, y)
#contr1 = plt.triplot(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1,\
       #    'k-',lw = 0.20)
#contr = plt.tripcolor(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1, strain[:, 0],\
         # edgecolors='k',cmap='rainbow' )
#contr = plt.tripcolor(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1, Dispx, shading='gouraud')
#contr = ax.triplot(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1, 'r-', lw = 1 )
contr = plt.tricontourf(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1, Dispy, 40, cmap='coolwarm' )
plt.colorbar()
print('The code took %.2f seconds to complete the run'%(time.time()-startTime) )
#ax.set_title('Contour plot of displacement')
#ax.set_xlabel('X ')
#ax.set_ylabel('Y')
#plt.show()
#for ii in range(nele):
 #   xvec, yvec = np.array([]), np.array([])
  #  dxvec, dyvec = np.array([]), np.array([])
   # elem = ElemData[ii]
    #for jj in range(len(elem)):
     #   xvec = np.append(xvec, NodeData[elem[jj]-1, 1])
      #  yvec = np.append(yvec, NodeData[elem[jj]-1, 2])
       # dxvec = np.append(dxvec, DeformedCoord[elem[jj]-1, 0])
        #dyvec = np.append(dyvec, DeformedCoord[elem[jj]-1, 1])
    #xvec = np.append(xvec, xvec[0])
    #yvec = np.append(yvec, yvec[0])
    #dxvec = np.append(dxvec, dxvec[0])
    #dyvec = np.append(dyvec, dyvec[0])
    #plt.plot(xvec, yvec, lw = '1', color = 'blue', label = 'undeformed')
    #plt.plot(dxvec, dyvec, '--', lw = '1', color = 'red', label ='deformed')
    #plt.legend(loc=0)
    
    