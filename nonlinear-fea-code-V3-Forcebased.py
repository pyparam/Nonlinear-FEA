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
    elem = ElemData[elabel]-1
    xi, xj, xm = NodeData[elem[0],1], NodeData[elem[1],1], NodeData[elem[2],1]#, 0, 20, 0
    yi, yj, ym = NodeData[elem[0],2], NodeData[elem[1],2], NodeData[elem[2],2]
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
        #enum = 
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
def NlGlobalAssembly(ElemData, NodeData, LoadNode, E, nu, rextglobal, kglobal, delU):
    """Assembly of tri element stiffness (ke) and nodal force vector (rext) into 
    the global form rglobl, kglobal respectively.
    rext: elemental nodal force vector
    ke: elemental stiffness matrix 
    kglobal: global stiffness matrix 
    """
    for el in range(nele):
        elem = ElemData[el]
        #enum = 
        B, A = Bmatrix(el, ElemData, NodeData, centroids)
        D = TangentMatrix(E, nu)
        tempM = np.dot(np.transpose(B),D)
        ke = A*np.dot(tempM,B) 
        stressi = A*LinearComputeStressStrain(ElemData, nele, NodeData, delU, E, nu, el)[1]
        stressi = stressi.reshape(3,1)
        connect = np.repeat(elem, 2)
        r = np.array([2, 1, 2, 1, 2, 1])
        rext = list(LoadNode[elem[0]-1,:]) + list(LoadNode[elem[1]-1,:])+ list(LoadNode[elem[2]-1,:])
        rint = np.matmul(np.transpose(B), stressi)
        #print(rext)
        for i in range(0, elemdof):
            #rextglobal[(2*connect[i]-r[i]), 0] = rextglobal[(2*connect[i]-r[i]), 0] + rext[i]
            rintglobal[(2*connect[i]-r[i]), 0] = rintglobal[(2*connect[i]-r[i]), 0] + rint[i]
            rextglobal[(2*connect[i]-r[i]), 0] = rext[i]
            #print(rglobal)
            #print('\n')
            for j in range(0,elemdof):
                kglobal[(2*connect[i]-r[i]), (2*connect[j]-r[j])] = kglobal[(2*connect[i]-r[i]), (2*connect[j]-r[j])] + \
                ke[i, j]
                #print(kglobal[(2*connect[i]-r[i]), (2*connect[j]-r[j])])
                #print(rglobal)
                #print('\n')
    rglobal = rextglobal-rintglobal
    return kglobal, rglobal, rextglobal
def DeAssemblyGlobalMatrix(kglobal, rglobal, SupportBCnodes, xdisp):
    kreduced, rreduced = kglobal.copy(), rglobal.copy()
    print('Displacement boundary condition successfully applied!')
    rreduced = np.ma.asarray(rreduced)
    kreduced = np.ma.asarray(kreduced)
    for icount in SupportBCnodes:
        rreduced[icount, 0] = np.ma.masked
        kreduced[icount, icount] = np.ma.masked 
    kreduced = np.ma.mask_rowcols(kreduced)
    rreduced =rreduced.compressed()
    kreduced =kreduced.compressed()
    kreduced = kreduced.reshape(len(rreduced), len(rreduced))
    
    return kreduced, rreduced
def SolveU(kreduced, rreduced, UnknownBCnodes, SupportBCnodes, xdisp):
    unknownU = np.linalg.solve(kreduced,rreduced)
    return unknownU

def assembleU(DisplacewithBC,UnknownBCnodes, SupportBCnodes, unknownU, xdisp):
    Displacement = DisplacewithBC.copy()
    print('Solution completed for unknown displacements')
    for indx in range(len(UnknownBCnodes)):
        index = int(UnknownBCnodes[indx])
        Displacement[index] =unknownU[indx] 
    for indx in range(len(SupportBCnodes)):
        index = int(SupportBCnodes[indx])
        Displacement[index] = 0#unknownU[indx] 
    return Displacement
def LinearComputeStressStrain(ElemData, nele, NodeData, delU, E, nu, e):
    strain, stress = np.zeros([nele, 4]), np.zeros([nele, 4])
    #for e in range(len(ElemData)):
    enum = ElemData[e]-1
        # elemDisp = [enum[0], enum[1], enum[2]]   
    B, Area = Bmatrix(e, ElemData, NodeData, centroids)     
    uarray = np.transpose(np.array([delU[2*enum[0]], delU[2*enum[0]+1],\
                                    delU[2*enum[1]], delU[2*enum[1]+1],\
                                    delU[2*enum[2]], delU[2*enum[2]+1]]))
    strain = np.matmul(B, uarray)   
    D = TangentMatrix(E, nu)
    stress = np.matmul(D, strain)
    return strain, stress    




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
rextglobal = np.zeros([nnode*ndof, 1])
rintglobal = np.zeros([nnode*ndof, 1])
strain, stress = np.zeros([nele, 4]), np.zeros([nele, 4])
print('Total no. of nodes : ', nnode)
print('Total no. of elements : ', nele)
print('Applying force on the nodes')
print('Elementwise solution started: Building global force and stiffness matrix')
################################################################################
DisplacewithBC = np.zeros([ndof*nnode,1])
xdisp = 1 # prescribed displacement BC
for i in range(len(NodeCoord)):
        if NodeCoord[i, 0] == 20.0:
            LoadNode[i, 0] =200000
print('Applying displacement boundary condition to the constrained node')
whereBC = [20, 0]
def GetNodeset(whereBC, DisplacewithBC):
    SupportBCnodes = np.array([]).astype(int)
    UnknownBCnodes = np.array([]).astype(int)
    for i in NodeintLabel:
        if NodeCoord[i, 0] ==0.0:#NodeCoord[i, 0] <= 5.0 and NodeCoord[i, 1] == 0.0 : #and 
            DisplacewithBC[2*i, 0] = 0
            DisplacewithBC[2*i+1, 0] = 0
            SupportBCnodes = np.append(SupportBCnodes, [2*i, 2*i+1])
#else:
        UnknownBCnodes = np.append(UnknownBCnodes, [2*i, 2*i+1])
    UnknownBCnodes = np.setxor1d(UnknownBCnodes, SupportBCnodes)
    return SupportBCnodes, UnknownBCnodes, DisplacewithBC
SupportBCnodes, UnknownBCnodes, DisplacewithBC =  GetNodeset(whereBC, DisplacewithBC)  


for idisp in range(0, 1):
    #SupportBCnodes, LoadBCnodes, UnknownBCnodes, DisplacewithBC =  getNodeset(whereBC, DisplacewithBC)  
    #DispMat = InitalDisp(DisplacewithBC)
    du = DisplacewithBC.copy()#DispInc(DispMat, idisp)
    #print(du.shape)
    du = du.reshape(len(du))
    kglobal, rglobal, rext = NlGlobalAssembly(ElemData, NodeData, LoadNode, E, nu, rextglobal, kglobal, du)
    kreduced, rreduced = DeAssemblyGlobalMatrix(kglobal, rglobal, SupportBCnodes, xdisp)
    solvedU =SolveU(kreduced, rreduced, UnknownBCnodes, SupportBCnodes, xdisp)
    #solvedU = solvedU.reshape(len(solvedU), 1)
    tempsolvedU = np.zeros(len(solvedU))#.reshape(len(solvedU), 1)
    for iconv in range(45):
        convcriteria = (sum(rglobal**2))/(1+ sum(rextglobal**2))
        print(convcriteria)
        tempsolvedU = tempsolvedU + solvedU
        #print(du)
        du = assembleU(DisplacewithBC,UnknownBCnodes, SupportBCnodes, tempsolvedU, xdisp)
        du = du.reshape(len(du))
       # print(du)
        kglobal, rglobal, rext = NlGlobalAssembly(ElemData, NodeData, LoadNode, E, nu, rextglobal, kglobal, du)
        kreduced, rreduced = DeAssemblyGlobalMatrix(kglobal, rglobal, SupportBCnodes, xdisp)
        solvedU =SolveU(kreduced, rreduced, UnknownBCnodes, SupportBCnodes, xdisp)
        
   # print(du)
    #print('\n')
tempsolvedU = tempsolvedU + solvedU
Displacement = assembleU(DisplacewithBC,UnknownBCnodes, SupportBCnodes, tempsolvedU, xdisp)   









############################################################################################################
DeformedCoord = np.zeros([nnode, 2])
print('Obtaining the deformed shape')    
for ii in range(nnode):
    DeformedCoord[ii, :] = NodeCoord[ii] + np.array([Displacement[2*ii, 0], Displacement[2*ii+1, 0]])
force = np.matmul(kglobal, Displacement)    
DeformNodeData= np.zeros([len(DeformedCoord), 3])
DeformNodeData[:, 0] = np.linspace(1, len(DeformedCoord), len(DeformedCoord))-1
DeformNodeData[:, 1] = DeformedCoord[:, 0]
DeformNodeData[:, 2] = DeformedCoord[:, 1]  
Dispx, forceX = np.array([]), np.array([])
Dispy, forceY = np.array([]), np.array([])
for i in range(nnode):
    Dispx = np.append(Dispx, Displacement[2*i])
    #print(Displacement[2*i])
   # print(2*i)
    Dispy = np.append(Dispy, Displacement[2*i+1])
    forceX = np.append(forceX, force[2*i])
    forceY = np.append(forceY, force[2*i+1])

plt.tricontourf(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1, forceX, 40, cmap='coolwarm' )
contr1 = plt.triplot(DeformedCoord[:,0], DeformedCoord[:,1], ElemData-1,\
           'k-',lw = 0.20)
plt.triplot(NodeCoord[:,0], NodeCoord[:,1], ElemData-1,\
           'k-',lw = 0.30)
plt.colorbar()



################################################################################
