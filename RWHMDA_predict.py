# encoding: UTF-8
import xlrd
import cmath
import math
from math import e
import numpy as np
import numpy.linalg as LA
import time
import scipy as sp

f1 = open(r'knowndiseasemicrobeinteraction.txt', 'r')
f3 = open(r'diseasesimilarity.txt', 'r')
N = 292
M = 39
# microbe, N
# disease, M
# N * M

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

def JH_simM(KM):
    JHM = KM.copy()
    return JHM

def JH_simD(KD, Disease_SSmatrix):
    JHD = 0.5*(KD+Disease_SSmatrix)
    return JHD
def weighted_incident_matrix(integrated_matrix_M,adjcent_Matrix):
    temp_M = adjcent_Matrix.copy()
    temp_MiRNA_diag = sum(integrated_matrix_M)
    Final_M = np.transpose(temp_M) *temp_MiRNA_diag
    Final_M = np.transpose(Final_M)
    return Final_M
def Prediction_Microbe(Interaction,Transtion):
    Predictionmatrix=np.zeros((N,M))
    for i in range(M):
        error = 1000000
        der = sum(Interaction[:,i])
        if der == 0:
            der += 1
        else:
            der = der
        initial = Interaction[:,i] / der

        initial0 = initial.copy()
        count=0
        while error > 1e-6:
            temp_initial=0.8*np.dot(initial,Transtion)+0.2*initial0
            error = LA.norm(temp_initial-initial,1)

            initial = temp_initial
            count+=1

        Predictionmatrix[:,i] = initial
    return Predictionmatrix

def for_transition_input(adjcent_M):
    Gauss_D_temp = Gauss_D(adjcent_M)
    Gauss_M_temp = Gauss_M(adjcent_M)
    Integrated_D_temp = JH_simD(Gauss_D_temp, F3)
    Integrated_M_temp = JH_simM(Gauss_M_temp)
    weighted_hyperedge_matrix_We_temp = np.diag(np.array([1./M]*M))
    node_Degree_matrix_Dv_temp = np.diag(np.dot(adjcent_M, sum(weighted_hyperedge_matrix_We_temp)))
    incident_matrix_W_temp = weighted_incident_matrix(Integrated_M_temp, adjcent_M)
    hyperedge_Degree_matrix_Dve_temp = np.diag(sum(incident_matrix_W_temp))
    Dv_inverse_temp = np.zeros((N, N))
    for i in range(N):
        Dv_inverse_temp[i][i] = 1. / (node_Degree_matrix_Dv_temp[i][i])
    Dve_inverse_temp = np.zeros((M, M))
    for j in range(M):
        if hyperedge_Degree_matrix_Dve_temp[j][j] ==0:
            Dve_inverse_temp[j][j] = 0
        else:
            Dve_inverse_temp[j][j] = 1. / (hyperedge_Degree_matrix_Dve_temp[j][j])
    P_part1_temp = np.dot(Dv_inverse_temp, adjcent_M)
    P_part2_temp = np.dot(weighted_hyperedge_matrix_We_temp, Dve_inverse_temp)
    temp_P_part3_here = np.dot(P_part1_temp, P_part2_temp)
    transition_matrix_P_here = np.dot(temp_P_part3_here, np.transpose(incident_matrix_W_temp))
    return transition_matrix_P_here
F1 = np.zeros((N, M))
for line in f1:
    arr = map(int, line.split())
    arr = list(arr)
    if arr[1] != 0:
        F1[arr[0] - 1][arr[1] - 1] = 1
print ('F1=', F1)                                           #F1 
print (F1.shape)
print (np.sum(F1))
F3 = np.zeros((M,M))

for line in f3:
    arr = map(float, line.split())
    arr = list(arr)
    if arr[1] != 0:
        F3[int(arr[0]) - 1][int(arr[1]) - 1] = arr[2]
                                                           #F3 
print (F3.shape)


Trans_Matrix = for_transition_input(F1)
Pre = Prediction_Microbe(F1,np.transpose(Trans_Matrix))
print(Pre)
'''
data =xlrd.open_workbook(r'disease-numbers.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
disease_data=[]
for row_index in range(table.nrows):
    disease_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8').decode('utf-8')])
print (disease_data)

data =xlrd.open_workbook(r'microbe-numbers.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
microbe_data=[]
for row_index in range(table.nrows):
    microbe_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8').decode('utf-8')])
print (microbe_data)
fw = open('prediction.txt','w')
str3=''
for i in range(N):
    for j in range(M):
        if F1[i][j]==0:
            str3+=str(disease_data[j][1])+'\t'+str(microbe_data[i][1])+'\t'+str(Pre[i][j].real)+'\n'
fw.write(str3)
fw.close()

print ('ok')


# np.savetxt('pre180108.txt',Pre)
'''