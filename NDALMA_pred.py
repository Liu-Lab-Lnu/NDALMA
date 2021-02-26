#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: UTF-8

import xlwt
import xlrd
import math
from math import e
import numpy as np
import time
import random
import operator
import pandas as pd
from sklearn.metrics import roc_auc_score

#Read similarity matrix
M = 770; N = 275;K_CONECT = 4966
Y = np.zeros((M ,N))
S_r = np.zeros((M ,M))
S_p = np.zeros((N ,N))
rFunctionalArray=np.zeros((M ,M))
pSemanticArray=S_p = np.zeros((N ,N))

f_MD1 = xlrd.open_workbook('./data/final_mi_seq_similarity_matrix.xlsx')
f_MD2 = xlrd.open_workbook('./data/final_lnc_seq_similarity_matrix.xlsx')

fw_PredictResult= open("./predictedresult.txt","w")


for i in range(M):
    for j in range(M):
        rFunctionalArray[i][j]=float(f_MD2.sheet_by_name('final_lnc_seq_similarity_matrix').cell_value(i,j))
for i in range(N):
    for j in range(N):
        pSemanticArray[i][j]=float(f_MD1.sheet_by_name('final_mi_seq_similarity_matrix').cell_value(i,j)) 


#Read interaction matrix
f_MDI = xlrd.open_workbook('./data/interaction_pair_seq.xlsx')
for i in range(K_CONECT):
    Y[int(f_MDI.sheet_by_name('one_zero_matrix').cell_value(i,0))][int(f_MDI.sheet_by_name('one_zero_matrix').cell_value(i,1))] =1

#Gaussian interaction profile kernel similarity for lncRNAs
Gama_d1 = 1. #Set Gama_d1 to 1.

GuassD=np.zeros((N,N))
Gama_d=N*Gama_d1/sum(sum(Y * Y)) 
for i in range(N):
            for j in range(N):
                GuassD[i][j]=e**(-Gama_d*(np.dot(Y[:,i]-Y[:,j],Y[:,i]-Y[:,j])))
#Gaussian interaction profile kernel similarity for miRNAs
Gama_m1 = 1.
GuassR=np.zeros((M,M))
Gama_m=M*Gama_m1/sum(sum(Y * Y)) 
for i in range(M):
            for j in range(M):
                GuassR[i][j]=e**(-Gama_m*(np.dot(Y[i,:]-Y[j,:],Y[i,:]-Y[j,:])))
#End of Gaussian similarity

#Integration of lncRNA similarity, miRNA similarity with Gaussian similarity
#lncRNA
for i in range(N):
    for j in range(N):
        if pSemanticArray[i][j] == 0:
            S_p[i][j] = GuassD[i][j]
        else:
            S_p[i][j] = pSemanticArray[i][j]
#miRNA
for i in range(M):
    for j in range(M):
        if rFunctionalArray[i][j] == 0:
            S_r[i][j] = GuassR[i][j]
        else:
            S_r[i][j] = rFunctionalArray[i][j]
                          
#lncRNA
#Step 1 
#Network distance computation and adjustment
D2=np.zeros((M,M))
D1=np.zeros((M,M))
D1=1/S_r
for i in range(len(S_r)):
    for j in range(len(S_r)):
        D2[i][j]=D1[i][j]/(np.sqrt(D1.sum(1)[i]*D1.sum(1)[j]))
        
#Step 2  
list1=np.zeros((N))
w=np.zeros((M,N))
U1=np.zeros((M,N))
for p in range(N):
    list1=Y[:,p]
    for i in range(len(D2)):
        if len(D2[i][list1==1]):
            w[i][p]=np.mean(D2[i,:])-np.mean(D2[i][list1==1])
        else:
            w[i][p]=np.mean(D2[i,:])-D2[i,0]

#Step 3  
#Calculate Score 
arr1=Y.sum(0)
arr2=Y.sum(1)
D3=np.zeros((M,N))
for j in range(N):
    for i in range(M):
        D3[i][j]=arr1[j]/(np.where(w[:,j]>=w[i][j],1,0).sum())

#miRNA
#Step 1  
#Network distance computation and adjustment
D4=np.zeros((N,N))
D5=np.zeros((M,N))
D6=np.zeros((N,N))
D6=1/S_p        
for i in range(len(S_p)):
    for j in range(len(S_p)):
        D4[i][j]=D6[i][j]/(np.sqrt(D6.sum(1)[i]*D6.sum(1)[j]))
#Step 2  
list2=np.zeros((M))        
r=np.zeros((M,N))
for p in range(M):
    list2=Y[p,:]
    for j in range(N):
        if len(D4[j][list2==1]):
            U1[p][j]=np.mean(D4[0:M,j])-np.mean(D4[j][list2==1])
        else:
            U1[p][j]=np.mean(D4[0:M,j])

#Step 3  
#Calculate Score 
for i in range(M):
    for j in range(N):
        D5[i][j]=arr2[i]/(np.where(U1[i,:]>=U1[i][j],1,0).sum())        
        
#Average Score 
P=np.zeros((M,N))
P=(D3+D5)/2  
##Save prediction results (column 1: miRNA name, column 2: lncRNA name, column 3: score)
Exc = []
for i in range(M):
    for j in range(N):
        if Y[i][j] != 1:
            ListNew = [P[i][j],j,i]
            Exc.append(ListNew)
Exc.sort(reverse=True)
for i in range(len(Exc)):
    fw_PredictResult.writelines([str(Exc[i][1]),"\t",str(Exc[i][2]),"\t",str(Exc[i][0]),"\n"])
fw_PredictResult.close()



