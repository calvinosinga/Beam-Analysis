import Gaussians as gs
import numpy as np
import matplotlib.pyplot as plt

fi = np.loadtxt('D:/data_CasAs_20171017/bl_108/day_0_i.txt')
fr = np.loadtxt('D:/data_CasAs_20171017/bl_108/day_0_r.txt')
dim = fi.shape
seconds = np.zeros(dim[0])
vis = gs.visibilities(fi,fr)
temp = np.zeros([dim[1],dim[0]])
for i in range(dim[0]):
    for j in range(dim[1]):
        temp[j][i]=vis[i][j]
print(temp.shape)
for s in range(len(seconds)):
    seconds[s]=s*60

degrees = gs.seconds_to_degrees(seconds, 0,temp)
print(degrees)

