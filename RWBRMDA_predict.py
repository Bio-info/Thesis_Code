# encoding: UTF-8
import xlrd
import cmath
import math
import random
from math import e
import numpy as np
import time
from numpy import linalg as LA
import scipy as sp
from scipy import optimize
from scipy.optimize import fmin
from scipy.linalg.misc import norm
from scipy.optimize import minimize
f1 = open(r'knowndiseasemirnainteraction.txt', 'r')
f2 = open(r'mirnasimilarity.txt', 'r')
f21 = open(r'weight_mirna.txt', 'r')
f3 = open(r'diseasesimilarity1.txt', 'r')
f31 = open(r'weight_disease.txt', 'r')
f4 = open(r'diseasesimilarity2.txt', 'r')
N = 495
M = 383
K = 15
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

def JH_simM(adj_matrix, RNA_FSmatrix, RNA_FSjqmatrix):

    KM = Gauss_M(adj_matrix)
    JHM = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            if RNA_FSjqmatrix[i][j] == 1:
                JHM[i][j] = RNA_FSmatrix[i][j]
            elif RNA_FSjqmatrix[i][j] == 0:
                JHM[i][j] = KM[i][j]
    return JHM

def JH_simD(adj_matrix, Disease_SSmatrix, Disease_SSjqmatrix):
    KD = Gauss_D(adj_matrix)
    JHD = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            if Disease_SSjqmatrix[i][j] == 1:
                JHD[i][j] = Disease_SSmatrix[i][j]
            elif Disease_SSjqmatrix[i][j] == 0:
                JHD[i][j] = KD[i][j]
    return JHD

def Lap_D(D_jhsim,adj_Matrix):
    Identity_D=np.identity(M)
    mD=sum(np.transpose(D_jhsim))
    mD=mD**(-0.5)
    D_dis=np.diag(mD)
    Laplas_D=Identity_D-np.dot(np.dot(D_dis,D_jhsim),D_dis)
    XD_Lap=np.dot(adj_Matrix,Laplas_D)
    return XD_Lap

def Lap_M(M_jhsim,adj_Matrix):
    Identity_M = np.identity(N)
    mM=sum(np.transpose(M_jhsim))
    mM=mM**(-0.5)
    D_mir = np.diag(mM)
    Laplas_M = Identity_M - np.dot(np.dot(D_mir , M_jhsim), D_mir )
    XM_Lap = np.dot(Laplas_M,adj_Matrix)
    return XM_Lap

F1 = np.zeros((N, M))
for line in f1:
    arr = map(int, line.split())
    if arr[1] != 0.:
        F1[arr[0] - 1][arr[1] - 1] = 1.
print 'F1=', F1                                           #F1 邻接矩阵
print F1.shape


F = []
lines = f2.readlines()
for j in range(len(lines)):
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F2 = np.array(F)                                           #F2  mRNA 功能类似矩阵
print 'F2=',F2
F = []
lines = f21.readlines()
for j in range(len(lines)):
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F21 = np.array(F)                                           #F21 mRNA 加权矩阵
print 'F21=',F21
print sum(F21)
F = []
lines = f3.readlines()
for j in range(M):
    if j % 1 == 0:
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F3 = np.array(F)                                            #F3 疾病语义类似性1

F = []
lines = f31.readlines()
for j in range(M):
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F31 = np.array(F)                                           #F31  疾病加权矩阵

F = []
lines = f4.readlines()
for j in range(M):
        line = [float(x) for x in lines[j].split()]
        F.append(line)
F4 = np.array(F)                                            #F4 疾病语义类似性2

F5=(1./2)*(F3+F4)                                           #F5 疾病语义类似性（两个相加除以2）


jh_D=JH_simD(F1,F5,F31)
print jh_D
print jh_D==np.transpose(jh_D)
jh_M=JH_simM(F1,F2,F21)
print jh_M
print jh_M == np.transpose(jh_M)

# 构造miRNA网络上的随机游走

def Prediction_Mirna (miRNA_func_matrix,Transtion):
    Predictionmatrix = np.zeros((N,N))

    for i in range(N):
        error = 1000000
        der = sum(miRNA_func_matrix[:,i])
        initial = miRNA_func_matrix[:,i] / der
        initial0 = initial.copy()
        count=0
        while error > 1e-6:
            temp_initial=0.9*np.dot(Transtion,initial)+0.1*initial0
            error = LA.norm(temp_initial-initial,1)
            initial = temp_initial
            count+=1
        Predictionmatrix[:,i] = initial
    return Predictionmatrix
# 根据功能类似性构造转移矩阵
start1 = time.clock()
Tem_for_normalize = sum(jh_M)

for i in range(Tem_for_normalize.shape[0]):
        Tem_for_normalize[i] = 1./Tem_for_normalize[i]

Tran_Matrix = jh_M * Tem_for_normalize
print Tran_Matrix
print sum(Tran_Matrix)

Pre = Prediction_Mirna(jh_M,Tran_Matrix) #对每个RNA随机游走之后的得分矩阵
print Pre
# np.savetxt('Pre.txt',Pre)
Score_matrix = np.zeros((N,M))# 对每个疾病的每个RNA打分矩阵
for i in range(M):
    Eigenvector_Matrix = np.zeros((N, 3))
    for j in range(N):
        tem_dic1 = [] # 疾病关联RNA的信息
        tem_dic2 = []
        if F1[j][i] == 1:
            unite = tuple((j,F1[j][i]))
            tem_dic1.append(unite)
        temp_score_list = list(Pre[:,j])
        score_list = list(enumerate(temp_score_list))
        score_list.sort(key = lambda x: x[1], reverse = True)
        First_k_list = score_list[0:K]
        First_k_listcopy = First_k_list[:]
        sum_1_weight = 0
        sum_all_weight = 0
        for x in tem_dic1:
            for y in First_k_listcopy:
                if x[0] == y[0]:
                    sum_1_weight += y[1]
        for z in First_k_listcopy:
            sum_all_weight += z[1]
        sum_0_weight = sum_all_weight - sum_1_weight
        Eigenvector = np.array([1,sum_1_weight,sum_0_weight])
        Eigenvector_Matrix[j] = Eigenvector

    def func(derta):
        sum = 0
        for p in range(N):
            sum += ( F1[p][i] * np.dot(derta,Eigenvector_Matrix[p]) - math.log(1 + e**np.dot(derta,Eigenvector_Matrix[p])))
        return 1./sum
    initial = []
    for q in range(3):
        initial.append(random.random())
    initial = np.array(initial)
    # opt_w = optimize.fmin_cg(func, initial)
    res = minimize(func,initial,method='BFGS')
    opt_w = res.x
    for r in range(N):
        prediction_score = math.exp(np.dot(opt_w,Eigenvector_Matrix[r]))/float(1+math.exp(np.dot(opt_w,Eigenvector_Matrix[r])))
        Score_matrix[r,i] = prediction_score
end1 = time.clock()
print Score_matrix
print 'running time is %s'%(end1- start1)
        # set1 = set(tem_dic1)
        # set2 = set(First_k_listcopy)
        # temp_set = set1 & set2
        # number_1 = len(temp_set)
#np.savetxt('Score_matrix.txt',Score_matrix)


data =xlrd.open_workbook(r'diseasenumbers.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
disease_data=[]
for row_index in range(table.nrows):
    disease_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])
print disease_data

data =xlrd.open_workbook(r'miRNAnumbers.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)

miRNA_data=[]
for row_index in range(table.nrows):
    miRNA_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])
print miRNA_data
f7w=open(r'results_for_all.txt','w')
str3=''
for i in range(N):
    for j in range(M):
        if F1[i][j]==0:
            str3+=str(disease_data[j][1])+'\t'+str(miRNA_data[i][1])+'\t'+str(Score_matrix[i][j].real)+'\n'
f7w.write(str3)
f7w.close()

print 'ok'


#data =xlrd.open_workbook(r'results_for_all.xlsx')
#table=data.sheets()[0]
#x=table.nrows
#for row_index in range(table.nrows):
    #f=open(r'%s.txt'%table.row(row_index)[0].value.encode('utf-8'),'a')
    #f.write(str(table.row(row_index)[0].value.encode('utf-8'))+'\t'+str(table.row(row_index)[1].value.encode('utf-8'))+'\t'+str(table.row(row_index)[2].value)+'\n')
#print 'okky'



