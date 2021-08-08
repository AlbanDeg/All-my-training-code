import matplotlib.pyplot as plt
import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set_theme()
import matplotlib.animation as animation
import sys
from matplotlib.backend_bases import MouseButton
##initialisation and input of all parameters

s='The program uses Python 3.9.6'
print(s.center(50))
print('')
s1='All these library are required to run '
print(s1.center(50))
print('')
s2='numpy,matplotlib,seaborn'
print(s2.center(50))

#map size
x_axis_length=52
y_axis_length=52

#library check
lib=['numpy','matplotlib.pyplot','seaborn']
for i in range(len(lib)):
    if lib[i] in sys.modules:
        pass
    else:
        print(' Could not run because a library required is missing :',lib[i])
        sys.exit()
print(' ')
print('library check : ok ')
print('loading...')
print('')
print('map size set on '+ str( x_axis_length-2)+'x'+str( y_axis_length-2))
print('WARNING : Huge amount of generation can be long to run so keep the amount above 1000')
print('')
s='Enter the number of generation '
print(s.center(50))
print('')
nb_gen=input()
s='Choose starting pattern : blinker ,block ,boat,random'
print(s.center(50))
print('')
starting_pattern=input()

#empty map
empty_map=np.zeros(shape=(y_axis_length,x_axis_length))

# #plot first empty map for test

# plt.figure(figsize=[10,10])
# ax = sns.heatmap(empty_map,cbar=False,yticklabels=False,xticklabels=False)
# plt.show()


#starting pos
first_gen_xpos=25
first_gen_ypos=25


##creating starting pattern

##oscillators
#blinker
if starting_pattern=='blinker':
    empty_map[first_gen_xpos,first_gen_ypos]=1
    empty_map[first_gen_xpos+1][first_gen_ypos+1]=0
    empty_map[first_gen_xpos-1][first_gen_ypos-1]=0
    empty_map[first_gen_xpos+1][first_gen_ypos-1]=0
    empty_map[first_gen_xpos-1][first_gen_ypos+1]=0
    empty_map[first_gen_xpos][first_gen_ypos+1]=1
    empty_map[first_gen_xpos][first_gen_ypos-1]=1
    empty_map[first_gen_xpos-1][first_gen_ypos]=0
    empty_map[first_gen_xpos+1][first_gen_ypos]=0
else:
    pass

#beacon
if starting_pattern=='beacon':
    empty_map[first_gen_xpos,first_gen_ypos]=1
    empty_map[first_gen_xpos+1][first_gen_ypos+1]=0
    empty_map[first_gen_xpos-1][first_gen_ypos-1]=0
    empty_map[first_gen_xpos+1][first_gen_ypos-1]=1
    empty_map[first_gen_xpos-1][first_gen_ypos+1]=1
    empty_map[first_gen_xpos][first_gen_ypos+1]=0
    empty_map[first_gen_xpos][first_gen_ypos-1]=1
    empty_map[first_gen_xpos-1][first_gen_ypos]=1
    empty_map[first_gen_xpos+1][first_gen_ypos]=1
else:
    pass

##stable patern
#blcok
if starting_pattern=='block':
    empty_map[first_gen_xpos,first_gen_ypos]=1
    empty_map[first_gen_xpos+1][first_gen_ypos+1]=0
    empty_map[first_gen_xpos-1][first_gen_ypos-1]=0
    empty_map[first_gen_xpos+1][first_gen_ypos-1]=1
    empty_map[first_gen_xpos-1][first_gen_ypos+1]=0
    empty_map[first_gen_xpos][first_gen_ypos+1]=0
    empty_map[first_gen_xpos][first_gen_ypos-1]=1
    empty_map[first_gen_xpos-1][first_gen_ypos]=0
    empty_map[first_gen_xpos+1][first_gen_ypos]=1
else:
    pass

#boat
if starting_pattern=='boat':
    empty_map[first_gen_xpos,first_gen_ypos]=0
    empty_map[first_gen_xpos+1][first_gen_ypos+1]=0
    empty_map[first_gen_xpos-1][first_gen_ypos-1]=0
    empty_map[first_gen_xpos+1][first_gen_ypos-1]=1
    empty_map[first_gen_xpos-1][first_gen_ypos+1]=1
    empty_map[first_gen_xpos][first_gen_ypos+1]=1
    empty_map[first_gen_xpos][first_gen_ypos-1]=0
    empty_map[first_gen_xpos-1][first_gen_ypos]=1
    empty_map[first_gen_xpos+1][first_gen_ypos]=1
else:
    pass
##random pattern
if starting_pattern=='random':
    empty_map=np.random.randint(2,size=(y_axis_length,x_axis_length))

#runing several generation

current_map=empty_map
all_map=[]
n=0
while n<int(nb_gen):
    new_map=np.zeros(shape=(y_axis_length,x_axis_length))
    for i in range(1,51):
        for j in range(1,51):
            env=np.sum([current_map[i+1][j+1],current_map[i+1][j],current_map[i][j+1],current_map[i-1][j-1],current_map[i][j-1],current_map[i-1][j],current_map[i-1][j+1],current_map[i+1][j-1]])
            if env==3 and current_map[i][j]==0:
                new_map[i][j]=1
            elif env==2 and current_map[i][j]==1 :
                new_map[i][j]=1
            elif env==3 and current_map[i][j]==1 :
                new_map[i][j]=1
            elif env==2 and current_map[i][j]==0 :
                new_map[i][j]=0
            elif env==3 and current_map[i][j]==0 :
                new_map[i][j]=1
            elif env<2 or env>3:
                new_map[i][j]=0
            else:
                pass
    current_map=new_map
    all_map.append(new_map)
    n=n+1

fig=plt.figure(figsize=[7,7])
def init():

    sns.heatmap(data=np.zeros(shape=(y_axis_length,x_axis_length)),cbar=False,yticklabels=False,xticklabels=False)
def animate(i):
    current_data=all_map[i]
    plt.title("Game of Life v2.0 generation n°"+str(i+1))
    sns.heatmap(data=current_data,cbar=False,yticklabels=False,xticklabels=False,cmap="YlGnBu")
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(all_map), repeat = False)
plt.show()

