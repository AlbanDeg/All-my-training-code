
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import seaborn as sns; sns.set_theme()
from matplotlib import cm


#first order sinusoidale perturbation with approximation

#paramter

w_fi=10
t=10
A=1/100
#plotting
w=np.linspace(-15,15,1000000)
p_fi=A*((np.sin((w_fi+w)*t/2))**2/((w_fi+w)/2)**2 + (np.sin((w_fi-w)*t/2))**2/((w_fi-w)/2)**2)

plt.figure()
plt.title("Perturbation sinosoïdale au premier ordre (approx) ")
plt.plot(w,p_fi,'r')
plt.xlabel('Pulsation de la perturbation w')
plt.ylabel('Probabilité de transition f-->i')
plt.show()

