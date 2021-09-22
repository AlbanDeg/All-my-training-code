import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import seaborn as sns; sns.set_theme()
from matplotlib import cm
from matplotlib.ticker import LinearLocator

one_D='on'
two_D='on'
#################################################################################################
## 1D
#paramater
if one_D=='on':
    y=[]
    x=[]
    t=[]
    L=10
    dom=np.linspace(-L/2,L/2,300)
    for point in dom:
        x.append(point)
    Norm=math.sqrt(2/L)
    n=12

    #ploting
    plt.figure()
    plt.xlabel('Position')
    plt.ylabel('Well')
    plt.title(' InfQuantum well PT : Wave function and probability of presence')
    plt.legend()
    def phi_impair(x):
        return Norm*math.cos(n*math.pi/L*x)
    def phisquared_impair(x):
        return (Norm*math.cos(n*math.pi/L*x))**2
    def phi_pair(x):
        return Norm*math.sin(n*math.pi/L*x)
    def phisquared_pair(x):
        return (Norm*math.sin(n*math.pi/L*x))**2
    if n%2==0:
        for i in x:
            y.append(phi_pair(i))
            t.append(phisquared_pair(i))
        plt.plot(x,t,color='black',label='Probability')
        plt.plot(x,y,'r',label='Wave function')
        plt.legend(loc="upper right")
        plt.show()
    elif n%2!=0:
        for i in x:
            y.append(phi_impair(i))
            t.append(phisquared_impair(i))
        plt.plot(x,t,color='black',label='Probability')
        plt.plot(x,y,'g',label='Wave function')
        plt.legend(loc="upper right")
        plt.show()

##############################################################################################################
## 2D
if two_D=='on':
    L_x=10
    L_y=10

    n_x=2
    n_y=6
    x=np.arange(-L_x/2,L_x/2,0.1)
    y=np.arange(-L_y/2,L_y/2,0.1)
    z_grid=np.zeros(shape=(100,100))
    def prob(x,y):
        return (math.sqrt(2/L_x)*math.sin(n_x*math.pi*x/L_x)*math.sqrt(2/L_y)*math.sin(n_y*math.pi*y/L_y))**2
    for i in range(len(x)):
        for j in range(len(y)):
            z_grid[i][j]=prob(x[i],y[j])
    plt.figure()
    plt.xlabel('X postion')
    plt.ylabel('Y position')
    plt.title(f'2D Quantum Well : Porbability of presence nx={n_x} ; ny={n_y}')
    sns.heatmap(data=z_grid,yticklabels=False,xticklabels=False,cmap="bwr")
    plt.show()


