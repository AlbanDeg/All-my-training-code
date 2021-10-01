import math
import numpy as np

    #config elec dans une seule couche
print("Entrer la configuration electronique, exemple:'d2' ")
config=str(input())

#toutes les couches elec format : (lettre, l , nb electrons dans la couche)
couches=[('s',0,2),('p',1,6),('d',2,10),('f',3,14),('g',4,18)]

# nb electrons dans le cas etudié
number_elec=int(config[1])
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
electron_settled=0
i=len(current_couche)-1
while True:
    if current_couche[i][1]==0:
        current_couche[i][1]='+'
        electron_settled=electron_settled+1

    if electron_settled==number_elec:
        break
    if current_couche[i][2]==0:
        electron_settled=electron_settled+1
        current_couche[i][2]='-'

    if electron_settled==number_elec:
        break
    i=i-1
M_L=0
#calcul de M_L
for j in range(len(current_couche)):
    for k in range(len(current_couche[j])):
        if current_couche[j][k]=='+' or current_couche[j][k]=='-':
            M_L=M_L+ current_couche[j][0]
#Calcul de M_S
M_S=number_elec*1/2

#Affichage
print(f"Pour la configuration {config} : ")
print(f"-{M_L} <= M_L <= {M_L}")
print(f"-{M_S} <= M_S <= {M_S}")


#Création du tableau originale
tab_ml_ms={}
tab_ms={}
all_ms=np.arange(-M_S,M_S,1)
all_ml=np.arange(-M_L,M_L,1)
for i in range(len(all_ml)):
        tab_ml_ms[f"ml={str(all_ml[i])}"]=tab_ms
for j in range(len(all_ms)):
    tab_ms[f"ms={str(all_ms[j])}"]={}
