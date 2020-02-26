#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 20:52:32 2019

@author: paramveer
"""
import numpy as np
NodeData = np.array([])    
ElemDatawithElemNo =   np.array([])  
#f2 = open('Elem-11.txt', 'w')
#f1 = open('Node-11.txt', 'w')
with open('Job-11.inp') as file:
    for line in file:
        if line[0:5] =='*Node':
           #2 print(line)
            for line in file:
                if line[0:8] !='*Element':
                    ndata = np.fromstring(line, dtype = np.float, sep = ',')
                    NodeData = np.append(NodeData, ndata)
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
                    ElemDatawithElemNo = np.append(ElemDatawithElemNo, edata)
                   # f2.write(line)
                else:
                    break
            #f2.close()
        
        
        
def ReadAbaqusInp():
    NodeDataInp = np.array([])    
    ElemDatawithElemNoInp =   np.array([])   
    with open('Job-1.inp') as file:
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
                                    #f2.close()
    return NodeDataInp, ElemDatawithElemNoInp           
         
print('Reading geometry and mesh data!')
NodeData, ElemDatawithElemNo = ReadAbaqusInp()