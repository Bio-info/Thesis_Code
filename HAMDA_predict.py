# encoding: GBK
import xlrd
import math
from math import e
import numpy as np
import time


f1=open(r'G:\Python程序\交叉验证\原始\1.数据\1.miRNA-disease关联数据\knowndiseasemirnainteraction.txt','r')
f2=open(r'G:\Python程序\交叉验证\原始\1.数据\4.miNA功能类似性矩阵\miRNA功能类似性矩阵.txt','r')
f21=open(r'G:\Python程序\交叉验证\原始\1.数据\4.miNA功能类似性矩阵\miRNA功能类似性加权矩阵.txt','r')
f3=open(r'G:\Python程序\交叉验证\原始\1.数据\2.disease semantic similarity 1\疾病语义类似性矩阵1.txt','r')
f31=open(r'G:\Python程序\交叉验证\原始\1.数据\2.disease semantic similarity 1\疾病语义类似性加权矩阵1.txt','r')
f4=open(r'G:\Python程序\交叉验证\原始\1.数据\3.disease semantic similarity 2\疾病语义类似性矩阵2.txt','r')
f41=open(r'G:\Python程序\交叉验证\原始\1.数据\3.disease semantic similarity 2\疾病语义类似性加权矩阵2.txt','r')
f0=open(r'G:\Python程序\交叉验证\RNA_degree.txt','w')
f00=open(r'G:\Python程序\交叉验证\Disease_degree.txt','w')
f1w=open(r'G:\Python程序\交叉验证\RS邻接矩阵.txt','w')
f2w=open(r'G:\Python程序\交叉验证\RR类似性矩阵.txt','w')
f3w=open(r'G:\Python程序\交叉验证\DD类似性矩阵.txt','w')
f4w=open(r'G:\Python程序\交叉验证\修改后的RNA类似性矩阵.txt','w')
f5w=open(r'G:\Python程序\交叉验证\Final_Result.txt','w')
f6w=open(r'G:\Python程序\交叉验证\S1.txt','w')
f7w=open(r'G:\Python程序\交叉验证\综合结果.txt','w')
#行为miRNA, N 
#列为disease, M 
#N为行数  M为列数

N=495
M=383
F1=np.zeros((495,383))

#先处理miRNA与disease的邻接矩阵
RNA_degree=[0 for i in range (N)]
Disease_degree=[0 for j in range(M)]

for line in f1:
    arr=map(int,line.split())# 把原始数据存在数组里
    #print arr
    if arr[1]!=0.:
        F1[arr[0]-1][arr[1]-1]=1.
    if arr[1]!=0.:
        RNA_degree[arr[0]-1]+=1
    if arr[0]!=0.:
        Disease_degree[arr[1]-1]+=1
#f0.write(str(RNA_degree)+'\t')
#f00.write(str(Disease_degree)+'\t')
print 'F1=',F1
print RNA_degree
print Disease_degree

#for i in range(N):
    #for j in range(M):
        #f1w.write(str(F1[i][j])+'\t')
    #f1w.write('\n')
#f1w.close()  #把邻接矩阵保存下来

#再处理miRNA之间的类似性矩阵
F=[]

lines=f2.readlines()
for j in range(N):
    if j%1==0:
        line=[float(x) for x in lines[j].split()]
        F.append(line)
F2= np.array(F) 
#把列表转换成可以做运算的np.ndarray形式
#print F2
#for i in range(N):
    #for j in range(N):
        #f2w.write(str(F2[i][j])+'\t')
    #f2w.write('\n')
#f2w.close()  # 把数据保存下来
F=[]
lines=f21.readlines()
for j in range(N):
    if j%1==0:
        line=[float(x) for x in lines[j].split()]
        F.append(line)
F21= np.array(F)

#再处理disease之间的类似性矩阵1
F=[]
lines=f3.readlines()
for j in range(M):
    if j%1==0:
        line=[float(x) for x in lines[j].split()]
        F.append(line)
F3= np.array(F)
#把列表转换成可以做运算的np.ndarray形式
#print F3
#for i in range(M):
    #for j in range(M):
        #f3w.write(str(F3[i][j])+'\t')
    #f3w.write('\n')
#f3w.close()  # 把数据保存下来
F=[]
lines=f31.readlines()
for j in range(M):
    if j%1==0:
        line=[float(x) for x in lines[j].split()]
        F.append(line)
F31= np.array(F)

def JH_simM(adj_matrix,RNA_FSmatrix,RNA_FSjqmatrix):
    KM=np.zeros((N,N))    
    JHM=np.zeros((N,N))        
    rm=N*1./sum(sum(adj_matrix*adj_matrix))
    for i in range(N):
        for j in range(N):
            KM[i][j]=e**(-rm*(np.dot(adj_matrix[i,:]-adj_matrix[j,:],adj_matrix[i,:]-adj_matrix[j,:])))
            if RNA_FSjqmatrix[i][j]==1:
                JHM[i][j]=RNA_FSmatrix[i][j]
            elif RNA_FSjqmatrix[i][j]==0:
                JHM[i][j]=KM[i][j]                
    return JHM
    
def JH_simD(adj_matrix,Disease_SSmatrix,Disease_SSjqmatrix):
    KD=np.zeros((M,M))
    JHD=np.zeros((M,M))
    T=adj_matrix.transpose()
    rd=M*1./sum(sum(T*T))
    for i in range(M):
        for j in range(M):
            KD[i][j]=e**(-rd*(np.dot(T[i]-T[j],T[i]-T[j])))
            if Disease_SSjqmatrix[i][j]==1:
                JHD[i][j]=Disease_SSmatrix[i][j]
            elif Disease_SSjqmatrix[i][j]==0:
                JHD[i][j]=KD[i][j]
    return JHD
JHD=JH_simD(F1,F3,F31)       
JHM=JH_simM(F1,F2,F21)
print JHD
print JHM
                
F4=np.zeros((N,N))
G=np.zeros((N,N))
H=np.zeros((N,N))
F32=np.zeros((N,N))
F5=np.zeros((M,M))
F6=np.zeros((1,M))
S=np.zeros((N,N))

start=time.clock()
for i in range(N):    
    for j in range(N):
        F32=JHD 
        H[i][j]=sum(F1[i,:])*sum(F1[j,:])
        F5=F32*F1[j,:]
        F6=np.dot(F1[i,:],F5)
        G[i][j]=sum(F6)
        if H[i][j]==0:
            S[i][j]=0
        else:
            S[i][j]=G[i][j]/H[i][j]        
 
end=time.clock()
print S
print 'Running time: %s Seconds'%(end-start)


alpha=0.7
lamda=0.8

S1=alpha * JHM + (1-alpha) * S 
print 'S1=',S1


P=np.zeros((N,N))
Q=np.zeros((N,N))
W=np.zeros((N,N))
D_degree=np.array(Disease_degree)
R_degree=np.array(RNA_degree)
Final_Result=np.zeros((N,M))
for i in range (N):
    for j in range(N):        
        Q[i][j]=(R_degree[i]**(1-lamda))*(R_degree[j]**lamda)        
        P[i][j]=S1[i][j]*sum(F1[i,:]*F1[j,:]/D_degree)
W=P/Q
Final_Result=np.dot(W,F1)
print 'Final_Result=',Final_Result
print 'OK'
#for i in range(N):
    #for j in range (M):
        
        #f5w.write(str(Final_Result[i][j])+'\t')
    #f5w.write('\n')
#f5w.close()


data =xlrd.open_workbook(r'G:\Python程序\交叉验证\原始\1.数据\1.miRNA-disease关联数据\disease数字编号.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
disease_data=[]
for row_index in range(table.nrows): 
    disease_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])    
print disease_data

data =xlrd.open_workbook(r'G:\Python程序\交叉验证\原始\1.数据\1.miRNA-disease关联数据\miRNA数字编号.xlsx')
table = data.sheets()[0]
print (table.nrows)
print (table.name)
x=table.nrows
miRNA_data=[]
for row_index in range(table.nrows): 
    miRNA_data.append([row_index+1,table.row(row_index)[1].value.encode('utf-8')])    
print miRNA_data


#for i in range(N):
    #for j in range(M):
        #if F1[i][j]==0:
            #f7w.write(str(disease_data[j][1])+'\t'+str(miRNA_data[i][1])+'\t'+str(Final_Result[i][j])+'\n')
            
#print 'ok'            
            
#data =xlrd.open_workbook(r'G:\Python程序\交叉验证\预测结果\所有疾病预测结果.xlsx')
#table=data.sheets()[0]
#x=table.nrows
#for row_index in range(table.nrows):
    #f=open(r'G:\Python程序\交叉验证\预测结果\%s.txt'%table.row(row_index)[0].value.encode('utf-8'),'a')
    #f.write(str(table.row(row_index)[0].value.encode('utf-8'))+'\t'+str(table.row(row_index)[1].value.encode('utf-8'))+'\t'+str(table.row(row_index)[2].value)+'\n')    
#print 'okky'        


def loocv(Neighbor_Matrix,W_matrix):
    P=np.zeros((N,N))
    Q=np.zeros((N,N))
    W=np.zeros((N,N))
    D_degree=np.array(Disease_degree)
    R_degree=np.array(RNA_degree)
    Final_Result=np.zeros((N,M))    
    for i in range(N):
        for j in range(N):
            Q[i][j]=(R_degree[i]**(1-lamda))*(R_degree[j]**lamda)
            P[i][j]=W_matrix[i][j]*sum(Neighbor_Matrix[i,:]*Neighbor_Matrix[j,:]/D_degree)
    W=P/Q
    Final_Result=np.dot(W,Neighbor_Matrix)  
    return Final_Result


def change(Neighbor_matrix,RNA_similarity,Disease_similarity):
    G=np.zeros((N,N))
    H=np.zeros((N,N))
    F5=np.zeros((M,M))
    F6=np.zeros((1,M))
    S=np.zeros((N,N))
    S1=np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            H[i][j]=sum(Neighbor_matrix[i,:])*sum(Neighbor_matrix[j,:])
            F5=Disease_similarity*Neighbor_matrix[j,:]
            F6=np.dot(Neighbor_matrix[i,:],F5)
            G[i][j]=sum(F6)
            if H[i][j]==0:
                S[i][j]=0
            else:
                S[i][j]=G[i][j]/H[i][j] 
    alpha=0.7
    lamda=0.8   
    S1=alpha * RNA_similarity + (1-alpha) * S    
    return S1  

def Ave_index(item,List):
    F=[]
    for i in range(len(List)):
        if List[i][2]==item:
            F.append(i+1)
    return float(sum(F))/len(F)

    
        
T=np.zeros((N,M))
O=np.zeros((N,M))
for i in range(N):
    
    start1=time.clock()
    for j in range(M):
        if F1[i][j]==1:
            F1[i][j]=0
            
            JHD=JH_simD(F1,F3,F31)       
            JHM=JH_simM(F1,F2,F21)
            FS1=change(F1,JHM,JHD)
            O=loocv(F1,FS1)
            print O
            l1=[]
            dic1=[]
            for r in range(N):
                if F1[r][j]==0:
                    unit1=tuple((r,j,O[r][j]))
                    dic1.append(unit1)
            l1=sorted(dic1,key=lambda item:item[2],reverse=True)
            b1=Ave_index(O[i][j],l1)
            print 'b1=',b1
            d1=open(r'G:\Python程序\交叉验证\验证结果\local排序位置.txt','a+')
            d1.write(str(i+1)+'\t'+str(j+1)+'\t'+str(b1)+'\n')        
            
            l=[]
            dic=[]
            for z in range(N):                                                
                for k in range(M):
                    if F1[z][k]==0:
                        unit=tuple((z,k,O[z][k]))
                        dic.append(unit)
            l=sorted(dic,key=lambda item:item[2],reverse=True) 
            #print l
            b=Ave_index(O[i][j],l)
            print 'b=',b
            d2=open(r'G:\Python程序\交叉验证\验证结果\global排序位置.txt','a+')
            d2.write(str(i+1)+'\t'+str(j+1)+'\t'+str(b)+'\n')           
            end1=time.clock()
            print 'Running time: %s Seconds'%(end1-start1)
            F1[i][j]=1
            

    