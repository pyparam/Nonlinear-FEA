#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 02:55:12 2019

@author: paramveer
"""
import numpy as np 
import matplotlib.tri as tri
import matplotlib.pyplot as plt

def ReadAbaqusInp(filename):
    NodeDataInp = np.array([])    
    ElemDatawithElemNoInp =   np.array([])   
    with open(filename) as file:
        for line in file:
            if line[0:5] =='*Node':
                #2 print(line)
                for line in file:
                    if line[0:8] !='*Element':
                        ndata = np.fromstring(line, dtype = np.float, sep = ',')
                        NodeDataInp = np.append(NodeDataInp, ndata)
                        #print(ndata)
                    else:
                        break
                        # f1.close()
            if line[0:8] =='*Element':
                    # if line[0:5] =='*Node':
                    # print(line)
                for line in file:
                    if line[0:1] !='*':
                        edata = np.fromstring(line, dtype = np.float, sep = ',')
                        ElemDatawithElemNoInp = np.append(ElemDatawithElemNoInp, edata)
                                # f2.write(line)
                    else:
                        break
    ## Reshaping from 1D array to numpy matrix
    NodeDataInp = NodeDataInp.reshape(int(len(NodeDataInp)/3), 3) #np.loadtxt('Node-221.txt', delimiter = ',')
    ElemDatawithElemNoInp = ElemDatawithElemNoInp.reshape(int(len(ElemDatawithElemNoInp)/4), 4)    
                           #f2.close()
    return NodeDataInp, ElemDatawithElemNoInp     
def MeshComputationfromPy():
    xlen = 20
    ylen = 10
    seed = 1
    xpoint,ypoint = np.array([]), np.array([])
    x1 = np.linspace(0,xlen,int(xlen//seed))
    y1 = 0*np.ones(int(xlen//seed))
    xpoint = np.append(xpoint, x1)
    ypoint = np.append(ypoint, y1)
    x2 = xlen*np.ones(int((ylen-seed)//seed))
    y2 = np.linspace(seed, ylen, int((ylen-seed)//seed))
    xpoint = np.append(xpoint, x2)
    ypoint = np.append(ypoint, y2)
    x3 = np.linspace(xlen-seed, 0, int((xlen-seed)//seed))
    y3 = ylen*np.ones(int((xlen-seed)//seed))
    xpoint = np.append(xpoint, x3)
    ypoint = np.append(ypoint, y3)
    x4 = 0*np.ones(int((ylen-2*seed)//seed))
    y4 = np.linspace(ylen-seed, seed, int((ylen-2*seed)//seed))
    xpoint = np.append(xpoint, x4)
    ypoint = np.append(ypoint, y4)
    xpoint, ypoint = xpoint, ypoint

    angle = 10
    radii = 3
    noOfangles = np.linspace(0, 360, int(360//angle))
    xcircle = radii*np.cos(noOfangles*np.pi/180.0)
    ycircle =radii*np.sin(noOfangles*np.pi/180.0)
    xpoint = np.append(xpoint, xcircle+10)
    ypoint = np.append(ypoint, ycircle+5)



    xrand = np.random.uniform(0, 20, 20000)
    yrand = np.random.uniform(0, 10, 20000)
    for i in range(len(xrand)):
        pos = (xrand[i]- 10)**2 + (yrand[i] - 5)**2
        distfromnode = ((xrand[i]- xpoint)**2 + (yrand[i]- ypoint)**2)**(1/2.0)
        if pos>radii**2 and min(distfromnode)>=.7*seed:
            xpoint=np.append(xpoint, xrand[i])
            ypoint=np.append(ypoint, yrand[i])
            #plt.scatter(xpoint, ypoint)

    triang = tri.Triangulation(xpoint,ypoint)
    triang.set_mask(np.hypot((xpoint[triang.triangles]+11).mean(axis=1), (ypoint[triang.triangles]+6).mean(axis=1)) < radii)
    
    
    trian = np.array([])

    checkTri = np.hypot((xpoint[triang.triangles]-10).mean(axis=1), (ypoint[triang.triangles]-5).mean(axis=1)) 
    for i in range(len(triang.triangles)):
        if  checkTri[i] > radii:
            trian = np.append(trian, triang.triangles[i])
            # print(len(trian))       
    trian = trian.reshape([int(len(trian)/3), 3])        
            #triang.set_mask(np.hypot((xpoint[triang.triangles]-10).mean(axis=1), (ypoint[triang.triangles]-5).mean(axis=1)) < radii)
   # plt.triplot(triang.x, triang.y, trian)

   # plt.triplot(triang.x, triang.y, triang.triangles)
   # plt.triplot(hh)
    
   
    nodex = triang.x
    print(len(nodex))
    nodey = triang.y
    nodelabel = np.arange(1, len(nodex)+1, 1)
    nodexy = np.zeros([len(nodex), 3])
    nodexy[:,0], nodexy[:,1], nodexy[:,2] = nodelabel, nodex, nodey
    #Elementxy = np.zeros([len(triang.triangles), 4])
    Elementxy = np.zeros([len(trian), 4])
    #Elementlabel = np.arange(0, len(triang.triangles), 1)
    Elementlabel = np.arange(1, len(trian)+1, 1)
    Elementxy[:,0], Elementxy[:,1:4]= Elementlabel, trian
    return nodexy, Elementxy
    
#print(ReadAbaqusInp('Job-Multihole.inp'))
#NodeData, ElemDatawithElemNo = ReadAbaqusInp('Job-1.inp')
#NodeData, ElemDatawithElemNo = MeshComputationfromPy()
#ElemDatawithElemNo = np.loadtxt('Elem-221.txt', delimiter = ',')
#ElemData  = ElemDatawithElemNo[:,1:4]
#ElemData = ElemData.astype(int)
#NodeCoord = NodeData[:, 1:3]
#plt.triplot(NodeCoord[:,0], NodeCoord[:,1], ElemData-1)


def MeshPlot(NodeXYCoord, ConnectMat, lc, ltype):
   # NodeCoord = NodeData[:, 1:3]
   # DeformedCoord = np.zeros([nnode, 2])
   # for ii in range(nnode):
        #DeformedCoord[ii, :] = NodeCoord[ii] + np.array(Displacement[2*ii], Displacement[2*ii+1])
    fig, ax = plt.subplots()
    noelem, nonode = np.shape(ConnectMat)
    for ii in range(noelem):
        xvec, yvec = np.array([]), np.array([])
        #dxvec, dyvec = np.array([]), np.array([])
        elem = ConnectMat[ii]
        for jj in range(len(elem)):
            xvec = np.append(xvec, NodeXYCoord[elem[jj]-1, 0])
            yvec = np.append(yvec, NodeXYCoord[elem[jj]-1, 1])
         #   dxvec = np.append(dxvec, DeformedCoord[elem[jj]-1, 0])
         #   dyvec = np.append(dyvec, DeformedCoord[elem[jj]-1, 1])
        xvec = np.append(xvec, xvec[0])
        yvec = np.append(yvec, yvec[0])
          #  dxvec = np.append(dxvec, dxvec[0])
          #  dyvec = np.append(dyvec, dyvec[0])

        ax.plot(xvec, yvec, ltype, lw = '1', color = lc, label = 'undeformed')
          #plt.plot(dxvec, dyvec, '--', lw = '1', color = 'red', label ='deformed')
    #plt.legend(loc=0)
    return ax
#MeshPlot(NodeCoord, ElemData, 'red', '-') 




