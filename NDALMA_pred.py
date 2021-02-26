#!/usr/bin/python
# -*- coding: UTF-8 -*-
# encoding: GBK

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

#将疾病语义类似性1文件存放于矩阵中
M = 771; N = 276;K_CONECT = 4966
Y = np.zeros((M ,N))
S_r = np.zeros((M ,M))
S_p = np.zeros((N ,N))
rFunctionalArray=np.zeros((M ,M))
pSemanticArray=S_p = np.zeros((N ,N))

f_MD1 = xlrd.open_workbook('./data/inal_mi_seq_similarity_matrix.xlsx')
f_MD2 = xlrd.open_workbook('./data/inal_lnc_seq_similarity_matrix.xlsx')

fw_PredictResult= open("./predictedresult.txt","w")


#将疾病语义类似性1文件存放于矩阵中
for i in range(M):
    for j in range(M):
        rFunctionalArray[i][j]=float(f_MD2.sheet_by_name('final_lnc_seq_similarity_matrix').cell_value(i,j))
for i in range(N):
    for j in range(N):
        pSemanticArray[i][j]=float(f_MD1.sheet_by_name('final_mi_seq_similarity_matrix').cell_value(i,j)) 
#for i in range(M):
#    for j in range(N):
#        b.append(Y[i][j])
#for i in range(13377):
#        fw_Y.writelines([str(b[i]),"\t"])
#fw_Y.close()
#高斯类似性矩阵计算
f_MDI = xlrd.open_workbook('./data/one_zero_matrix.xlsx')
for i in range(K_CONECT):
    Y[int(f_MDI.sheet_by_name('one_zero_matrix').cell_value(i,0))][int(f_MDI.sheet_by_name('one_zero_matrix').cell_value(i,1))] =1
        #高斯类似性矩阵计算
Gama_d1 = 1.
        #疾病高斯相似性
GuassD=np.zeros((N,N))
Gama_d=N*Gama_d1/sum(sum(Y * Y)) #新的带宽参数Gama_d1设置为1.
for i in range(N):
            for j in range(N):
                GuassD[i][j]=e**(-Gama_d*(np.dot(Y[:,i]-Y[:,j],Y[:,i]-Y[:,j])))
Gama_m1 = 1.
        #RNA高斯相似性
GuassR=np.zeros((M,M))
Gama_m=M*Gama_m1/sum(sum(Y * Y)) #新的带宽参数Gama_d1设置为1.
for i in range(M):
            for j in range(M):
                GuassR[i][j]=e**(-Gama_m*(np.dot(Y[i,:]-Y[j,:],Y[i,:]-Y[j,:])))
#高斯相似性计算完毕
#语义相似矩阵与高斯相似矩阵加权
#疾病加权计算
for i in range(N):
    for j in range(N):
        if pSemanticArray[i][j] == 0:
            S_p[i][j] = GuassD[i][j]
        else:
            S_p[i][j] = pSemanticArray[i][j]
        #miRNA加权计算
for i in range(M):
    for j in range(M):
        if rFunctionalArray[i][j] == 0:
            S_r[i][j] = GuassR[i][j]
        else:
            S_r[i][j] = rFunctionalArray[i][j]
#        for i in range(N):
#            for j in range(N):
#                if pSemanticArray[i][j] == 0:
#                    S_p[i][j] = GuassD[i][j]
#                else:
#                    S_p[i][j] = (pSemanticArray[i][j]+GuassD[i][j])/2
#                #miRNA加权计算
#        for i in range(M):
#            for j in range(M):
#                if rFunctionalArray[i][j] == 0:
#                    S_r[i][j] = GuassR[i][j]
#                else:
#                    S_r[i][j] = (GuassR[i][j]+rFunctionalArray[i][j])/2
#        for i in range(M):
#            for j in range(M):
#                S_r[i][j] = (GuassR[i][j]+rFunctionalArray[i][j])/2
#        for i in range(N):
#            for j in range(N):
#                S_p[i][j] = (pSemanticArray[i][j]+GuassD[i][j])/2
#        for i in range(M):
#            for j in range(M):    
#                S_r[i][j] = rFunctionalArray[i][j]
#                if S_r[i][j]==0:
#                    S_r[i][j]=0.1
#        for i in range(N):
#            for j in range(N):    
#                S_p[i][j] = pSemanticArray[i][j]
#                if S_p[i][j]==0:
#                    S_p[i][j]=0.1
                          
#第一步 
#调整距离
D2=np.zeros((M,M))
D1=np.zeros((M,M))
D1=1/S_r
for i in range(len(S_r)):
    for j in range(len(S_r)):
        D2[i][j]=D1[i][j]/(np.sqrt(D1.sum(1)[i]*D1.sum(1)[j]))
#第二步 
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
#        dm=0
#        for p in range(N):
#            list1=Y[:,p]
#            for o in range(len(list1)):
#                if list1[o]==1:
#                    dm+=1
#            for i in range(len(D2)):
#                if dm>0:
#                    w[i][p]=np.mean(D2[i,:])-np.mean(D2[i,:dm])
#                else:
#                    w[i][p]=np.mean(D2[i,:])
#            dm=0
#第三步
#计算得分
arr1=Y.sum(0)
arr2=Y.sum(1)
D3=np.zeros((M,N))
for j in range(N):
    for i in range(M):
        D3[i][j]=arr1[j]/(np.where(w[:,j]>=w[i][j],1,0).sum())
#第一步 
#调整距离
D4=np.zeros((N,N))
D5=np.zeros((M,N))
D6=np.zeros((N,N))
D6=1/S_p        
for i in range(len(S_p)):
    for j in range(len(S_p)):
        D4[i][j]=D6[i][j]/(np.sqrt(D6.sum(1)[i]*D6.sum(1)[j]))
#第二步 
list2=np.zeros((M))        
r=np.zeros((M,N))
for p in range(M):
    list2=Y[p,:]
    for j in range(N):
        if len(D4[j][list2==1]):
            U1[p][j]=np.mean(D4[0:M,j])-np.mean(D4[j][list2==1])
        else:
            U1[p][j]=np.mean(D4[0:M,j])
#        md=0
#        for p in range(M):
#            list2=Y[p,:]
#            for p in range(len(list2)):
#                md+=1
#            for j in range(N):
#                if md>0:
#                    U1[p][j]=np.mean(D4[:,j])-np.mean(D4[:md,j])
#                else:
#                    U1[p][j]=np.mean(D4[:,j])-D4[0,j]
#            md=0
#第三步
#计算得分
for i in range(M):
    for j in range(N):
        D5[i][j]=arr2[i]/(np.where(U1[i,:]>=U1[i][j],1,0).sum())        
#平均得分
P=np.zeros((M,N))
P=(D3+D5)/2  
##给出不包含已知关联的excel从大到小排序（1列疾病名，2列RNA名，3列分数）   
Exc = []
for i in range(M):
    for j in range(N):
    # if Y[i][j] != 1:
        ListNew = [P[i][j],j,i]
        Exc.append(ListNew)
Exc.sort(reverse=True)
for i in range(len(Exc)):
    fw_PredictResult.writelines([str(Exc[i][1]),"\t",str(Exc[i][2]),"\t",str(Exc[i][0]),"\n"])
fw_PredictResult.close()



