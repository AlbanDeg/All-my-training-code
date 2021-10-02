import math
import numpy as np
import random as rdm

    #config elec dans une seule couche
print("Entrer la configuration electronique, exemple:'d2' ")
config=str(input())

#toutes les couches elec format : (lettre, l , nb electrons dans la couche)
couches=[('s',0,2),('p',1,6),('d',2,10),('f',3,14),('g',4,18)]

# nb electrons dans le cas etudié


number_elec = 0
for char in config:
    if char.isalpha() ==True:
        number_elec=int(config.replace(char,""))

letter=config[0]
#longueur de la couche etudié
for couche in couches:
    if letter==couche[0]:
        len_couche=couche[2]


#dégénérésence
g=int(math.factorial(len_couche)/(math.factorial(number_elec)*math.factorial(len_couche-number_elec)))

#construction de la couche
current_couche={}
for i in range (int(len_couche/2)):
    current_couche[i]=[0,0,0]

#affectation des valeurs de L dans chaque case
if letter=="s":
    current_couche[0][0]=0
if letter=="p":
    current_couche[0][0]=-1
    for i in range(2):
        current_couche[i+1][0]=current_couche[i][0]+1
if letter=="d":
    current_couche[0][0]=-2
    for i in range(4):
        current_couche[i+1][0]=current_couche[i][0]+1
if letter=="f":
    current_couche[0][0]=-3
    for i in range(6):
        current_couche[i+1][0]=current_couche[i][0]+1
if letter=="g":
    current_couche[0][0]=-4
    for i in range(8):
        current_couche[i+1][0]=current_couche[i][0]+1

#Ml et Ms : On remplie la couche en maximisant L et en tenant compte du principe de Pauli

#remplissage pour le calcul de M_L
electron_settled=0
i=len(current_couche)-1
while electron_settled<number_elec:
    if current_couche[i][1]==0:
        current_couche[i][1]='+'
        electron_settled=electron_settled+1

    if electron_settled == number_elec:
        break

    if current_couche[i][2]==0:
        electron_settled=electron_settled+1
        current_couche[i][2]='-'

    if electron_settled == number_elec:
        break

    i=i-1

M_L=0
#calcul de M_L
for j in range(len(current_couche)):
    for k in range(len(current_couche[j])):
        if current_couche[j][k]=='+' or current_couche[j][k]=='-':
            M_L=M_L+ current_couche[j][0]

#reinitialisation de la couche pour le calcul de M_S
current_couche={}
for i in range (int(len_couche/2)):
    current_couche[i]=[0,0,0]

if letter=="s":
    current_couche[0][0]=0
if letter=="p":
    current_couche[0][0]=-1
    for i in range(2):
        current_couche[i+1][0]=current_couche[i][0]+1
if letter=="d":
    current_couche[0][0]=-2
    for i in range(4):
        current_couche[i+1][0]=current_couche[i][0]+1
if letter=="f":
    current_couche[0][0]=-3
    for i in range(6):
        current_couche[i+1][0]=current_couche[i][0]+1
if letter=="g":
    current_couche[0][0]=-4
    for i in range(8):
        current_couche[i+1][0]=current_couche[i][0]+1


M_S=0
#remplissage pour le calcul de M_S
if number_elec-len(current_couche)<0:
    for i in range(number_elec):
        current_couche[i][1]='+'
elif number_elec-len(current_couche)==0:
    M_S=0
elif number_elec-len(current_couche)>0:
    down=number_elec-len(current_couche)
    for i in range(len(current_couche)):
        current_couche[i][1]='+'
    for i in range(down):
        current_couche[i][2]='-'

#Calcul de M_S
for j in range(len(current_couche)):
    for k in range(len(current_couche[j])):
        if current_couche[j][k]=='+':
            M_S=M_S+0.5
        if current_couche[j][k]=='-':
            M_S=M_S-0.5






#Création du tableau originale
tab_ml_ms={}
tab_ms={}
all_ms=np.arange(-M_S,M_S+1,1)
all_ml=np.arange(-M_L,M_L+1,1)
#création de tout les l
all_l=np.arange(-current_couche[len(current_couche)-1][0],current_couche[len(current_couche)-1][0]+1,1)

## TEST NEW TAB ML MS
tab_MLMS={}
key_list=[]
for ml in all_ml:
    for ms in all_ms:
        key_list.append([ml,ms])

for key in key_list:
    tab_MLMS[key[0],key[1]]=[]

#remplissage du tableau
couche_test=[]
for i in range(len_couche):
    couche_test.append(0)

#création de la config de départ
electron_settled=0
i=0
while electron_settled<number_elec:
    couche_test[i]=1
    electron_settled=electron_settled+1
    i=i+1

all_state=[]
nb_state=0
while nb_state<g:
    current_config=rdm.sample(couche_test,len_couche)
    if current_config not in all_state:
        all_state.append(current_config)
        nb_state=nb_state+1
#création des det de slatter
# 1 == occupé
#0 === vide

all_slatter_det=[]
i=0
for state in all_state:
    slatter_det=[]
    for electron_pos in range(len(state)):
        if state[electron_pos]==1 and electron_pos<int(len_couche/2):
                slatter_det.append([all_l[electron_pos],'+'])
        elif state[electron_pos]==1 and electron_pos>=int(len_couche/2):
                slatter_det.append([all_l[electron_pos-int(len_couche/2)],'-'])
        elif state[electron_pos]==0:
            pass
        else:
            pass
    all_slatter_det.append(slatter_det)

#remplissage du tableau
#calcul du ml et ms pour chaque det de slatter
print(tab_ml_ms)
nb_ket=0
for s_det in all_slatter_det:
    Ml=0
    Ms=0
    a=[]
    b=[]
    for electron in s_det:
        b.append(electron[0])
        if electron[1]=='+':
            a.append(0.5)
        elif electron[1]=='-':
            a.append(-0.5)
    Ml=np.sum(b)
    Ms=np.sum(a)
    nb_ket=nb_ket+1

    tab_MLMS[Ml,Ms].append(s_det)

def All_det(a):
    if a=='ok':
        for det in all_slatter_det:
            print(det)


print("")
#Affichage
print(f"Pour la configuration {config}, g(nombre d'états)={g} : ")
print(f"-{M_L} <= M_L <= {M_L}")
print(f"-{M_S} <= M_S <= {M_S}")
print(f"Pour consulter le Tableau ML MS de la configuration {config}, ecrire : TabMLMS[ML,MS]       exemple TabMLMS[-4,0] pour d2")
print(f"Pour savoir le nombre de determinant de slatter d'un couple ML MS : len(TabMLMS[ML][MS])")
print(f"Pour voir les {g} det de Slatter : All_det('ok')")








