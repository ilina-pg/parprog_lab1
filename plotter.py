import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

print('LOADING')
dataframe = pd.read_csv('output.csv', sep=',',dtype={'u': float, 'x':float, 't':float})
print('END')
x = dataframe['x']
t = dataframe['t']
u = dataframe['u']

MAXVALUE = 2147483640

def filter(arr, spl):
    ans = [0 for i in range(len(arr)//spl)]
    for i in range(len(ans)):
        ans[i] = arr[i] if abs(arr[i]) < 0.15e6 else 0.15e6#arr[i-1]#2*np.sign(arr[i*spl])
    return ans

spl = 1
t = [i/100 for i in range(int(min(t)*100), int(max(t)*100 + 1))]
x = [i/100 for i in range(int(min(x)*100), int(max(x)*100 + 1))]
u = np.array(filter(list(u), spl))

print('FILTERED')

matrix = np.zeros((int(len(t)), int(len(x))))

for i in range(len(t)):
    for j in range(len(x)):
        #matrix[i][j] = u[i + j*len(t)]
        matrix[i][j] = u[j + i*len(x)]

'''for i in range(50):
    for j in range(50):
        matrix[i+100][j+100] = -10000000'''
        
print('MATRIX DONE')

fig = plt.figure(figsize=(12, 12))                    
ax_3d = fig.add_subplot(projection='3d')
xgrid, ygrid = np.meshgrid(x, t)
zgrid = matrix
ax_3d.plot_surface(xgrid, ygrid, zgrid, cmap='plasma')

ax_3d.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax_3d.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax_3d.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

plt.xlabel('x')
plt.ylabel('t')
#plt.zlabel('u')

ax_3d.tick_params(axis='both', which='major', labelsize=14)
ax_3d.tick_params(axis='both', which='minor', labelsize=12)

fig.savefig('plot.png', format='png', transparent = True, dpi=200)