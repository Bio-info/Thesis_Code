# encoding: UTF-8
import xlrd
import cmath
import math
from math import e
import numpy as np
import numpy.linalg as LA

N = 495
M = 383


def Transition_Matrix(Relation_Matrix):
    row = np.shape(Relation_Matrix)[0]
    col = np.shape(Relation_Matrix)[1]
    Transition_P = np.zeros((row,col))
    value, vector = LA.eig(Relation_Matrix)
    value = list(value)
    value = [item.real for item in value]
    max_lambda = max(value)
    max_index = value.index(max(value))
    max_vector = vector[:, max_index]
    max_vector = [u.real for u in list(max_vector)]
    temp = sum(max_vector)
    norm_max_vector = [x ** 0.5 / temp ** 0.5 for x in max_vector]
    for i in range(row):
        for j in range(col):
            if norm_max_vector[i] == 0.:
                Transition_P[i][j] = 0
            else:
                Transition_P[i][j] = (Relation_Matrix[i][j] * norm_max_vector[j]) / (max_lambda * norm_max_vector[i])
        Transition_P[i] = Transition_P[i]/(np.sum(Transition_P[i]))
    return Transition_P
def Prediction_MiRNA(Interaction,Transtion):
    row = np.shape(Interaction)[0]
    col = np.shape(Interaction)[1]
    Predictionmatrix=np.zeros((row,col))
    Final_Prediction = np.zeros((N,M))
    for k in range(row):
        error = 1000000
        der = sum(Interaction[k])
        if der == 0:
            der += 1
        else:
            der = der
        initial = Interaction[:,k] / der
        count=0
        while error > 1e-6:
            temp_initial=np.dot(np.transpose(Transtion),initial)
            error = LA.norm(temp_initial-initial,1)
            initial = temp_initial
            count+=1
        Predictionmatrix[k] = initial
        P_1 = Predictionmatrix[0:383,383:col]
        P_2 = Predictionmatrix[383:row,0:383]
        Final_Prediction = (np.transpose(P_1) + P_2) * 0.5
    return Final_Prediction

Adj_Matrix = np.zeros((N,M))     # 邻接矩阵
with open('knowndiseasemirnainteraction.txt') as File1:
    for line in File1:
        arr = map(int, line.split())
        arr = list(arr)
        if arr[1] != 0:
            Adj_Matrix[arr[0] - 1][arr[1] - 1] = 1
Mirna_Matrix = []            # RNA功能联系矩阵
with open('weight_mirna.txt') as File2:
    lines = File2.readlines()
    for j in range(N):
        line = [int(x) for x in lines[j].split()]
        Mirna_Matrix.append(line)
    Mirna_Matrix = np.array(Mirna_Matrix)
    for k in range(N):
        Mirna_Matrix[k][k] = 0.
Disease_Matrix = []           # 疾病语义联系矩阵
with open('weight_disease.txt') as File3:
    lines = File3.readlines()
    for j in range(M):
        line = [int(x) for x in lines[j].split()]
        Disease_Matrix.append(line)
    Disease_Matrix = np.array(Disease_Matrix)
    for k in range(M):
        Disease_Matrix[k][k] = 0.
temp1 = np.hstack((Disease_Matrix, np.transpose(Adj_Matrix)))
temp2 = np.hstack((Adj_Matrix, Mirna_Matrix))
Hyper_Matrix = np.vstack((temp1,temp2))    # 异质网络邻接矩阵


Transition_P_matrix = Transition_Matrix(Hyper_Matrix)
print(Transition_P_matrix)
print(np.sum(Transition_P_matrix,axis =1))
#np.savetxt('Transition_matrix.txt',Transition_P_matrix)
Prediction_as_MiRNA = Prediction_MiRNA(Hyper_Matrix,Transition_P_matrix)
print (Prediction_as_MiRNA)
print (np.shape(Prediction_as_MiRNA))
#np.savetxt('Prediction_matrix.txt',Prediction_as_MiRNA)
# print(sum(Mirna_Matrix))
# print(sum(Disease_Matrix))
# print(sum(Hyper_Matrix))
# print(Hyper_Matrix == np.transpose(Hyper_Matrix))
# value, vector = LA.eig(Hyper_Matrix)
# value = list(value)
# value = [item.real for item in value]
# print ('value = ', value)
# max_lambda = max(value)
# print ('max_lambda = ', max_lambda)
# max_index = value.index(max(value))
# max_vector = vector[:, max_index]
# print('before_abs=',max_vector)
# max_vector = [u.real for u in list(max_vector)]
# print('max_vector = ', max_vector)
# temp = sum(max_vector)
# print('sum_vector = ', temp)
# norm_max_vector = [x ** 0.5 / temp ** 0.5 for x in max_vector]
# print ('norm_max_vector = ',norm_max_vector)
# print (np.sum(np.dot(np.array(norm_max_vector),np.array(norm_max_vector))))
# print (np.sum(Mirna_Matrix - np.transpose(Mirna_Matrix)))