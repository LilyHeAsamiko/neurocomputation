# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 16:36:09 2020

@author: LilyHeAsamiko
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import preprocessing
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D

#regression on TF
dataset = pd.read_csv(r'/Users/he/Downloads/WiML-master/MAdata_DGinMCout.txt',sep='\t', header = 0)
origindata = np.exp(dataset)
dataset1 = pd.read_csv(r'/Users/he/Downloads/WiML-master/raw_data_DGinMCout.tsv',sep='\t', header = 0)
fadata = pd.read_csv(r'/Users/he/Downloads/WiML-master/Alignment_for_family_PTHR11679_SF35.txt',sep='\t', header = 0)

#X = origindata
#X.iloc[0:25,0] = dataset1.iloc[0:25,5]
#X.iloc[0:25,1] = dataset1.iloc[25:50,5]
#X.iloc[0:25,2] = dataset1.iloc[50:75,5]
X= np.zeros((np.shape(dataset)))
X[0:25,0] = dataset1.iloc[0:25,5]
X[0:25,1] = dataset1.iloc[25:50,5]
X[0:25,2] = dataset1.iloc[50:75,5]

N,c = np.shape(origindata)
#compare V1: dg_v, V2: dg_c
beta = np.dot(X.T,origindata)/N
beta = np.array(beta)
beta0 = beta
steps = 100
btemp = 0.1
db = 0.1
lse = np.ones(np.shape(origindata))*10**5
res = np.zeros((c,steps))
#regression
for b in range(c):
    for s in range(steps):
#        temp = btemp*X.iloc[:,b]
#        temp = np.unique((np.array(origindata.iloc[:,b]).reshape(N,1)-temp.reshape(N,1)-np.random.normal(0,1,N).reshape(N,1)))**2
#        temp = np.array(btemp*X[:,b])
#        temp = (np.array(origindata.iloc[:,b])-temp-np.random.normal(0,1,N).T)**2
        e = np.random.normal(0,1,N).T
        temp = (np.array(origindata.iloc[:,b])-np.array(btemp*X[:,b])-e)**2
        temp[np.isnan(temp)]=0
#        if sum(temp.reshape(N,1) - lse.reshape(N,1))<0:
        if sum(temp - lse[:,b])<0:
            beta[b,b] = btemp
            btemp += db
            lse[:,b] = temp 
            res[b,s] = np.sqrt(sum(lse[:,b])/N)
        else:
            btemp = (btemp +db)/2
            e = np.random.normal(0,1,N).T
            temp= (np.array(origindata.iloc[:,b])-np.array(btemp*X[:,b])-e)**2
            if sum(temp - lse[:,b])<0:
                beta[b,b] = btemp
                btemp += db
                lse[:,b] = temp 
                res[b,s] = np.sqrt(sum(lse[:,b])/N)
            else:
                res[b,s] = np.sqrt(res[b,s-1]/N)
         #       temp[np.isnan(temp)]=0
#        res[b,s] = np.mean(res[b,0：s-1]/s)
C = np.corrcoef(X,X)
plt.pcolor(C[0:25,0:25])
plt.title('correlation of three cells TF')
LD = np.zeros((np.shape(origindata)))
#LD score:
for i in range(c):
    LD[:,i] = origindata.iloc[:,i]**2+origindata.iloc[:,i]*origindata.iloc[:,np.mod(i+1,3)]+origindata.iloc[:,i]*origindata.iloc[:,np.mod(i+2,3)]
plt.pcolor(LD)
plt.title('LD_Score of three cells TF')
     
plt.scatter(origindata, np.dot(X,beta))
plt.title('three cells TF regression')


#Compare enrichment of fa dataset(as the sequence is not complete but only alignment data, there we compare enrichment only on family fasta sequence)
Str_Tnig = np.array(fadata.iloc[0:125],dtype = str)
Str_Onil = np.array(fadata.iloc[127:252],dtype = str)
#Common allele:
rows = np.shape(Str_Onil)[0]    
cols = len(str(Str_Onil[0,:]))
co = np.zeros((rows, cols))
cE = 0
S = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
for i in range(np.shape(Str_Onil)[0]):
    co[i,:] = str(Str_Onil[i,:])==str(Str_Tnig[i,:])
    cE += sum(co[i,:]>0)
M = rows*cols 
Enrichment_Tnig_Onil = np.zeros((rows, 2,2))
TnigID = np.zeros((rows, cols,len(S)))
OnilID = TnigID 
for i in range(rows):
    for s in range(len(S)):
#            TnigID[i,:,s] = str(Str_Tnig[i,:])[~co[i,:]]==S[s]
#        OnilID[i,:,s] = str(Str_Onil[i,~co[i,:]]==S[s]
#            OnilID[i,:,s] = str(Str_Onil[i,:])[~co[i,:]]==S[s]
        #temp1 = str(Str_Tnig[i,])[:]
        TnigID[i,:,s] = str(Str_Tnig[i,]).find(S[s])
        #print(len(temp1))
        #temp2 = str(Str_Onil[i,])[:]
        #print(len(temp2))
        OnilID[i,:,s] = str(Str_Onil[i,]).find(S[s])
        #tempTid = TnigID[i,:,s]>=0
        #tempOid = OnilID[i,:,s]>=0
        #tempc = np.corrcoef(TnigID[i,tempTid,s],OnilID[i,tempOid,s])
        #tempc[np.isnan(tempc)] = 0
        #Enrichment_Tnig_Onil[i,s,:,:] = tempc
    tempTid = TnigID[i,:,:]>=0
    tempOid = OnilID[i,:,:]>=0
    tempc = np.corrcoef(TnigID[i,tempTid],OnilID[i,tempOid])
    tempc[np.isnan(tempc)] = 0
    Enrichment_Tnig_Onil[i,:,:] = tempc

#fig = plt.figure()
#ax = fig.add_subplot(111,projection = '3d')
#ax.scatter(np.linspace(0,2,10),np.linspace(0,2,10),Enrichment_Tnig_Onil.T)
Enrichment1 = sum(Enrichment_Tnig_Onil[Enrichment_Tnig_Onil>0]**2)/cE/((sum(Enrichment_Tnig_Onil[Enrichment_Tnig_Onil<=0]**2))/(rows*cols-cE)+1)