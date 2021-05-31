# encoding: GBK
import xlrd
import cmath
import math
from math import e
import numpy as np
import time
import scipy as sp
from scipy import optimize
from scipy.optimize import fmin
from scipy.linalg.misc import norm
from scipy.optimize import minimize

f1=open(r'1.miRNA-disease关联数据\knowndiseasemirnainteraction.txt','r')
f2=open(r'4.miNA功能类似性矩阵\miRNA功能类似性矩阵.txt','r')
f3=open(r'2.disease semantic similarity 1\疾病语义类似性矩阵1.txt','r')
f4=open(r'3.disease semantic similarity 2\疾病语义类似性矩阵2.txt','r')
f7w=open(r'所有疾病预测结果.txt','w')
N=495
M=383
#行为miRNA, N 
#列为disease, M 
#N为行数  M为列数
def Gauss_M(adj_matrix):
    GM = np.zeros((N, N))
    rm = N * 1. / sum(sum(adj_matrix * adj_matrix))
    for i in range(N):
        for j in range(N):
            GM[i][j] = e ** (-rm * (np.dot(adj_matrix[i, :] - adj_matrix[j, :], adj_matrix[i, :] - adj_matrix[j, :])))
    return GM


def Gauss_D(adj_matrix):
    GD = np.zeros((M, M))
    T = adj_matrix.transpose()
    rd = M * 1. / sum(sum(T * T))
    for i in range(M):
        for j in range(M):
            GD[i][j] = e ** (-rd * (np.dot(T[i] - T[j], T[i] - T[j])))
    return GD


def AFinal_M(tM, tD, tQM, tQD):
    for i in range(len(tM)):
        tM[i] = tM[i].real
    for i in range(N):
        for j in range(N):
            tQM[i][j] = tQM[i][j].real
    C10 = np.kron(tD, tM)
    C20 = C10.copy()
    for i in range(len(C1)):
        C20[i] = 1. / (C10[i] + lamda * 1.)
    C30 = C20 * C10
    TtQM = tQM.transpose()
    TtQD = tQD.transpose()
    C40 = np.dot(TtQM, F1)
    C50 = np.dot(C40, tQD)
    C60 = []
    for i in range(M):
        for j in range(N):
            C60.append(C50[j][i])
    C60 = np.array(C60)
    C70 = C30.copy()
    for i in range(len(C70)):
        C70[i] = C70[i] * C60[i]
    C0 = np.zeros((N, M))
    for z in range(M):
        for k in range(N):
            C0[k][z] = C70[N * z + k]
    A00 = np.dot(tQM, C0)
    Anew = np.dot(A00, TtQD)
    return Anew



F1 = np.zeros((N, M))
for line in f1:
    arr = map(int, line.split())
    if arr[1] != 0.:
        F1[arr[0] - 1][arr[1] - 1] = 1.
print 'F1=', F1                                           #F1 邻接矩阵
print F1.shape
F = []
lines = f2.readlines()
for j in range(N):
    if j % 1 == 0:
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F2 = np.array(F)                                           #F2  mRNA 功能类似矩阵

F = []
lines = f3.readlines()
for j in range(M):
    if j % 1 == 0:
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F3 = np.array(F)                                            #F3 疾病语义类似性1

F = []
lines = f4.readlines()
for j in range(M):
    if j % 1 == 0:
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F4 = np.array(F)                                            #F4 疾病语义类似性2


GD = Gauss_D(F1)
GM = Gauss_M(F1)
KM0 = 1./2*(F2+GM)
KD0 = 1./3*(F3+F4+GD)
print 'KD0=', KD0
print 'KM0=', KM0
TKM0 = KM0.transpose()
TKD0 = KD0.transpose()
E = (KD0 + TKD0) / 2.
H = (KM0 + TKM0) / 2.
TD, QD = np.linalg.eig((KD0 + TKD0) / 2.)
TM, QM = np.linalg.eig((KM0 + TKM0) / 2.)


for i in range(len(TM)):
    TM[i] = TM[i].real
for i in range(N):
    for j in range(N):
        QM[i][j] = QM[i][j].real


lamda = 1.0
C1 = np.kron(TD, TM)
print 'C1=', C1
C2 = C1.copy()
for i in range(len(C1)):
    C2[i] = 1. / (C1[i] + lamda * 1.)
C3 = C2 * C1


TQM = QM.transpose()
TQD = QD.transpose()
C4 = np.dot(TQM, F1)
C5 = np.dot(C4, QD)
C6 = []
for i in range(M):
    for j in range(N):
        C6.append(C5[j][i])
C6 = np.array(C6)
print 'c6=', C6
C7 = C3.copy()
for i in range(len(C7)):
    C7[i] = C7[i] * C6[i]
print 'c3=', C3
print 'c7=', C7
print len(C7)
C = np.zeros((N, M))
for z in range(M):
    for k in range(N):
        C[k][z] = C7[N * z + k]
print C
print C.shape

A1 = np.dot(QM, C)
A = np.dot(A1, TQD)
print 'A=', A
for i in range(N):
    for j in range(M):
        A[i][j]=A[i][j].real
a = []
for i in range(M):
    for j in range(N):
        a.append(A[j][i].real)
a = np.array(a)
print 'a=', a
print len(a)



betaD0=[1./3,1./3,1./3]
betaD0=np.array(betaD0)
betaM0=[1./2,1./2]
betaM0=np.array(betaM0)
U=F1-A*(lamda/2.)

derta=0.5

def func0 (beta) :
    MD = beta[0] * np.dot(np.dot(KM0,A),F3) + beta[1] * np.dot(np.dot(KM0,A),F4) + beta[2] * np.dot(np.dot(KM0,A),GD)
    norm1=norm(U-MD,'fro')
    norm2 = norm(beta, 2)
    return 1./(2*lamda*N*M)*norm1+derta*norm2**2
initial=betaD0
cons=({'type':'eq','fun':lambda x: x[0]+x[1]+x[2]-1})
bnds=((0,1),(0,1),(0,1))
res=minimize(func0,initial,bounds=bnds,constraints=cons)
print res.x
print res.message
betaD0_new=res.x
print 'betaD0_new=',betaD0_new
KD1=betaD0_new[0]*F3+betaD0_new[1]*F4+betaD0_new[2]*GD
TKM0 = KM0.transpose()
TKD1 = KD1.transpose()
E1 = (KD1 + TKD1) / 2.
H1 = (KM0 + TKM0) / 2.
TD1, QD1 = np.linalg.eig(E1)
TM1, QM1 = np.linalg.eig(H1)
print QD1.shape
print QM1.shape
A_New=AFinal_M(TM1,TD1,QM1,QD1)
U_New=F1-A_New*(lamda/2.)
print 'A_New=',A_New

def func1(beta):
    MM = beta[0] * np.dot(F2, np.dot(A_New, KD1)) + beta[1] * np.dot(GM, np.dot(A_New, KD1))
    norm1=norm(U_New-MM,'fro')
    norm2=norm(beta,2)
    return 1./(2*lamda*N*M)*norm1+derta*norm2**2
initial=betaM0
cons=({'type':'eq','fun':lambda x: x[0]+x[1]-1})
bnds=((0,1),(0,1))
res=minimize(func1,initial,bounds=bnds,constraints=cons)
print res.x
print res.message
betaM0_new=res.x
print 'betaM0_new=',betaM0_new
KM1=betaM0_new[0]*F2+betaM0_new[1]*GM
TKM1 = KM1.transpose()
TKD1 = KD1.transpose()
E2 = (KD1 + TKD1) / 2.
H2 = (KM1 + TKM1) / 2.
TD2, QD2 = np.linalg.eig(E2)
TM2, QM2 = np.linalg.eig(H2)
A_New1=AFinal_M(TM2,TD2,QM2,QD2)
U_New1=F1-A_New1*(lamda/2.)
print 'A_New1=',A_New1

data =xlrd.open_workbook(r'1.miRNA-disease关联数据\disease数字编号.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
disease_data=[]
for row_index in range(table.nrows):
    disease_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])
print disease_data

data =xlrd.open_workbook(r'1.miRNA-disease关联数据\miRNA数字编号.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
miRNA_data=[]
for row_index in range(table.nrows):
    miRNA_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])
print miRNA_data

str3=''
for i in range(N):
    for j in range(M):
        if F1[i][j]==0:
            str3+=str(disease_data[j][1])+'\t'+str(miRNA_data[i][1])+'\t'+str(A_New1[i][j].real)+'\n'
f7w.write(str3)
f7w.close()

print 'ok'


#data =xlrd.open_workbook(r'所有疾病预测结果.xlsx')
#table=data.sheets()[0]
#x=table.nrows
#for row_index in range(table.nrows):
    #f=open(r'%s.txt'%table.row(row_index)[0].value.encode('utf-8'),'a')
    #f.write(str(table.row(row_index)[0].value.encode('utf-8'))+'\t'+str(table.row(row_index)[1].value.encode('utf-8'))+'\t'+str(table.row(row_index)[2].value)+'\n')    
#print 'okky' 