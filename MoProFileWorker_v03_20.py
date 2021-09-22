# importation des librairies
import numpy as np
import re
import os
from os import getcwd, chdir, mkdir
import matplotlib.pyplot as plt
import shutil
from operator import itemgetter
import math
from math import sqrt
import scipy
from scipy import stats
from sklearn.metrics import r2_score
import seaborn as sns
import sys
import time



# si on souhaite afficher(print dans la cmd) les verifications(afficher les bases et et la verification de l'orthogonalités) pour une configuration ab specifique : ---> mettre : verif_ab=1
verif_xy=0
verif_yx=0
verif_zx=0
verif_xz=0
verif_yz=0
verif_zy=0
verif_bxy=0
verif_byx=0
verif_bzx=0
verif_bxz=0
verif_byz=0
verif_bzy=0

s='The program uses Python 3.9.0'
print(s.center(50))
print('')
s1='All these library are required to run '
print(s1.center(50))
s2='numpy,os,re,matplotlib,shutil,operator,math,scipy,sklearn.metrics,seaborn'
print(s2.center(50))
find_type_atom_by_dip_qua_config='OFF'
print('Enter path to  your MoProViewer format (.par) data set :')
path_to_the_bigset=input()
print(' ')
print('Enter path to where output data files will be created !!! The atomic type file used must be in this path :')
path_to_txtfiles=input()
print(' ')
print('Enter path to an EMPTY DIRECTORY where your dataset with TPOL will be added')
desti=input()
print(' ')
print('Enter all graph/picture output path : ')
image_path=input()
print(' ')
print('First round with the program yes/no')
rep=input()
if rep=='yes':
    tpol='ON'
    conversion='ON'
    print(' ')



else:
    tpol='OFF'
    conversion='OFF'


print('atomic types file define : no or "filename.txt" : ')
rep=input()
if rep =='no':
    code_exit='OFF'
else:
    code_exit='ON'
    print(' ')
    print('path to empty ELMAM file ')
    elmam_path=input()
    path_to_liste=rep
start_time = time.time()

toogle_dim=0
create_distr_eig=1


##-------------------------Ouverture des fichiers------------------------------------------##
#library check
lib=['numpy','os','matplotlib.pyplot','shutil','operator','math','scipy','sklearn.metrics','seaborn','time']
for i in range(len(lib)):
    if lib[i] in sys.modules:
        pass
    else:
        print(' Could not run because a library required is missing :',lib[i])
        sys.exit()
print(' ')
print('library check : ok ')
print('loading...')
if conversion=='ON':
    os.chdir(path_to_the_bigset)
#repertoire ou est stocké le set [doit seulement contenir le set et aucun autre fichier!!]
    file = os.listdir()



    for x in range(len(file)):
        os.chdir(path_to_the_bigset)
        os.rename(file[x],file[x][:-5]+".txt")

#repertoire ou est stocké le set [doit seulement contenir le set et aucun autre fichier!!]
os.chdir(path_to_the_bigset)

fichiers = os.listdir()

for x in range(len(fichiers)):
    a=[]
    posatm=[]                                 #vecteur position des atomes
    ligne = []                                # création d'une liste de liste vide,
    os.chdir(path_to_the_bigset)
    with open (fichiers[x], "r") as f:
        for li in f :                         # pour toutes les lignes du fichier
            s = li.strip ("\n\r")             # on enlève les caractères de fin de ligne
            l = s.split (" ")                 # on découpe en colonnes
            ligne.append (l)
#nettoyage des chiffres


#fonction de nettoyage des lettres pour avoir les numeros des voisins
    def chiffre(X):
        return int(re.sub('[ABCDEFGHIJKLMNOPQRSTUVWXYZl_]' ,'',X))
#pour des cas speciaux ou le voisin est de la forme "2_C11" cette fonction enlève le "2_", cela ne fonctionne que si la condition d'unicité des numeros atomiques est respéctée dans tout les fichiers
    def residu(X):
        return X[2:]
##------------------Condition de marche du code sur les noms des atomes-------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if ligne[i][j]=="ATOM":
                if(int(ligne[i][1]) - chiffre(ligne[i][2])!= 0):
                    print('les noms ne sont pas uniques',fichiers[x])
                    sys.breakpoint()
                else:
                    pass

    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j]=="ATOM"):
    # si le mot ATOM est dans la ligne --> affectation des 3 coordonéees xyz dans un vecteur "posatm"
                posatm.append([ligne[i][7],ligne[i][8],ligne[i][9]])
    #le fichier texte ne renvoit que des type'str' il faut donc les changer en chiffre
    patm=np.array(posatm,dtype=float)
    #le vecteur patm comprend les 3 coordonées de tout les atomes du dossier
    o=[]
    v1=[]
    v2=[]
    #On introduit l'orthogonalisation qui retourne quelque soit la configuration géométrique les vecteurs dans une matrice (e_x,e_y,e_z) basé sur le procédé de Gram-Schmidt
    def ortho (atom,voisin1,voisin2,A):
        u1=[]
        u2=[]
        V1=[]
        V2=[]
        V1=np.array([voisin1[0]-atom[0],voisin1[1]-atom[1],voisin1[2]-atom[2]])
        V2=np.array([voisin2[0]-atom[0],voisin2[1]-atom[1],voisin2[2]-atom[2]])
        u2=np.array(V2-(np.vdot(V2,V1)/np.vdot(V1,V1))*V1) #Orthogonalisation de u2 selon u1
        e1=(V1/(np.linalg.norm(V1)))
        e2=((u2/(np.linalg.norm(u2))))
        e3=(np.cross(e1,e2)) # Produit vectoriel pour avoir le 3eme vecteur orthogonal

        #On place différentes conditions pour que l'algo redonne  les vecteurs dans le bon ordre
        if (A==12): #Configuration XY
            return np.array([e1,e2,e3])
        elif (A==13):  # XZ
            return np.array([e1,e3,e2])
        elif (A==21): # YX
            return np.array([e2,e1,e3])
        elif (A==23): #YZ
            return np.array([e3,e1,e2])
        elif (A==31): # ZX
            return np.array([e2,e3,e1])
        else: # ZY
            return np.array([e3,e2,e1])



##-----------------------------------------CONFIGURATION X Y --------------------------------------------##
    #trouver les repère dorientation XY
    o_xy=[]
    v1_xy=[]
    v2_xy=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("XY") ):
                o_xy.append(ligne[i-1][1])          #nom des atomes du fichier
                if(len(ligne[i][4])>3):
                    v1_xy.append(residu(ligne[i][4])) #nom des premiers voisins ex : 'C23'
                else:
                    v1_xy.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_xy.append(residu(ligne[i][8])) #nom des deuxièmes voisins ex : 'C24'
                else:
                    v2_xy.append(ligne[i][8])
    v_1v_xy=[]
    v_2v_xy=[]
    v_o_xy=[]
    for i in range(len(v1_xy)):      #on enelve les caracetres pour avoir juste les numeros des voisins ex 'C23' --> 23
        v_1v_xy.append(chiffre(v1_xy[i]))
    for i in range(len(v2_xy)):
        v_2v_xy.append(chiffre(v2_xy[i]))     #numero deuxieme vosin
    for i in range(len(o_xy)):
        v_o_xy.append(chiffre(o_xy[i]))       #numero atome etudié

    #Calcul et affichages des bases XY du dossier
    #attention!! les coordonées de latome numero N sont stockees dans le vecteur patm[N-1]
    bases_xy=[]
    for k in range(len(o_xy)):
        if (verif_xy==1):
            print("---------------------CONFIGURATION XY-------------------------")
            print("-----------------------Base atome",v_o_xy[k],"-----------------")
            print("                                                   ")
            print(ortho(patm[v_o_xy[k]-1],patm[v_1v_xy[k]-1],patm[v_2v_xy[k]-1],12))
            print("                                                   ")
        bases_xy.append([(ortho(patm[v_o_xy[k]-1],patm[v_1v_xy[k]-1],patm[v_2v_xy[k]-1],12))[0],(ortho(patm[v_o_xy[k]-1],patm[v_1v_xy[k]-1],patm[v_2v_xy[k]-1],12))[1],(ortho(patm[v_o_xy[k]-1],patm[v_1v_xy[k]-1],patm[v_2v_xy[k]-1],12))[2]])
    #boucle de verification de lorthogonalité des bases pour une certaine précision sur le determinant
    precision=0.00000001
    if (verif_xy==1):
        for t in range(len(bases_xy)):
            epsilon=np.linalg.det(bases_xy[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(abs(epsilon-1))
                print("la base de latome", v_o_xy[t] ," est orthonormee")


##------------------------------------------CONFIGURATION ZX------------------------------------------------------##
    o_zx=[]
    v1_zx=[]
    v2_zx=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("ZX") ):
                o_zx.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_zx.append(residu(ligne[i][4]))
                else:
                    v1_zx.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_zx.append(residu(ligne[i][8]))
                else:
                    v2_zx.append(ligne[i][8])
    v_1v_zx=[]
    v_2v_zx=[]
    v_o_zx=[]
    for i in range(len(v1_zx)):
        v_1v_zx.append(chiffre(v1_zx[i]))
    for i in range(len(v2_zx)):
        v_2v_zx.append(chiffre(v2_zx[i]))
    for i in range(len(o_zx)):
        v_o_zx.append(chiffre(o_zx[i]))

    bases_zx=[]
    for k in range(len(o_zx)):
        if (verif_zx==1):
            print("---------------------CONFIGURATION ZX-------------------------")
            print("------------------------Base atome",v_o_zx[k],"----------------")
            print("                                                   ")
            print(ortho(patm[v_o_zx[k]-1],patm[v_1v_zx[k]-1],patm[v_2v_zx[k]-1],31))
            print("                                                   ")
        bases_zx.append([ortho(patm[v_o_zx[k]-1],patm[v_1v_zx[k]-1],patm[v_2v_zx[k]-1],31)[0],ortho(patm[v_o_zx[k]-1],patm[v_1v_zx[k]-1],patm[v_2v_zx[k]-1],31)[1],ortho(patm[v_o_zx[k]-1],patm[v_1v_zx[k]-1],patm[v_2v_zx[k]-1],31)[2]])

    epsilon_t_zx=0

    precision=1e-15
    if (verif_zx==1):
        for t in range(len(bases_zx)):
            epsilon=np.linalg.det(bases_zx[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(abs(epsilon-1))
                print("la base de latome", v_o_zx[t] ," est orthonormee")

##-------------------------------------CONFIGURATION ZY----------------------------------##
    o_zy=[]
    v1_zy=[]
    v2_zy=[]
    A=32
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("ZY") ):
                o_zy.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_zy.append(residu(ligne[i][4]))
                else:
                    v1_zy.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_zy.append(residu(ligne[i][8]))
                else:
                    v2_zy.append(ligne[i][8])
    v_1v_zy=[]
    v_2v_zy=[]
    v_o_zy=[]
    for i in range(len(v1_zy)):
        v_1v_zy.append(chiffre(v1_zy[i]))
    for i in range(len(v2_zy)):
        v_2v_zy.append(chiffre(v2_zy[i]))
    for i in range(len(o_zy)):
        v_o_zy.append(chiffre(o_zy[i]))

    bases_zy=[]
    for k in range(len(o_zy)):
        if (verif_zy==1):
            print("---------------------CONFIGURATION ZY-------------------------")
            print("-----------------------Base atome",v_o_zy[k],"-----------------")
            print("                                                   ")
            print(ortho(patm[v_o_zy[k]-1],patm[v_1v_zy[k]-1],patm[v_2v_zy[k]-1],32))
            print("                                                   ")
        bases_zy.append([(ortho(patm[v_o_zy[k]-1],patm[v_1v_zy[k]-1],patm[v_2v_zy[k]-1],32))[0],(ortho(patm[v_o_zy[k]-1],patm[v_1v_zy[k]-1],patm[v_2v_zy[k]-1],32))[1],(ortho(patm[v_o_zy[k]-1],patm[v_1v_zy[k]-1],patm[v_2v_zy[k]-1],32))[2]])

    epsilon_t_zy=0
    precision=1e-15
    if (verif_zy==1):
        for t in range(len(bases_zy)):
            epsilon=np.linalg.det(bases_zy[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(abs(epsilon-1))
                print("la base de latome", v_o_zy[t] ," est orthonormee")


##------------------------------CONFIGURATION YZ --------------------------------------------##

    o_yz=[]
    v1_yz=[]
    v2_yz=[]
    A=32
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("YZ") ):
                o_yz.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_yz.append(residu(ligne[i][4]))
                else:
                    v1_yz.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_yz.append(residu(ligne[i][8]))
                else:
                    v2_yz.append(ligne[i][8])
    v_1v_yz=[]
    v_2v_yz=[]
    v_o_yz=[]
    for i in range(len(v1_yz)):
        v_1v_yz.append(chiffre(v1_yz[i]))
    for i in range(len(v2_yz)):
        v_2v_yz.append(chiffre(v2_yz[i]))
    for i in range(len(o_yz)):
        v_o_yz.append(chiffre(o_yz[i]))

    bases_yz=[]
    for k in range(len(o_yz)):
        if (verif_yz==1):
            print("---------------------CONFIGURATION YZ-------------------------")
            print("------------------------Base atome",v_o_yz[k],"----------------")
            print("                                                   ")
            print(ortho(patm[v_o_yz[k]-1],patm[v_1v_yz[k]-1],patm[v_2v_yz[k]-1],23))
            print("                                                   ")
        bases_yz.append([(ortho(patm[v_o_yz[k]-1],patm[v_1v_yz[k]-1],patm[v_2v_yz[k]-1],23))[0],(ortho(patm[v_o_yz[k]-1],patm[v_1v_yz[k]-1],patm[v_2v_yz[k]-1],23))[1],(ortho(patm[v_o_yz[k]-1],patm[v_1v_yz[k]-1],patm[v_2v_yz[k]-1],23))[2]])

    epsilon_t_yz=0

    precision=1e-15
    if (verif_yz==1):
        for t in range(len(bases_yz)):
            epsilon=np.linalg.det(bases_yz[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(abs(epsilon-1))
                print("la base de latome", v_o_yz[t] ," est orthonormee")


##---------------------------------------CONFIGURATION XZ -------------------------------##
    o_xz=[]
    v1_xz=[]
    v2_xz=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("XZ") ):
                o_xz.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_xz.append(residu(ligne[i][4]))
                else:
                    v1_xz.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_xz.append(residu(ligne[i][8]))
                else:
                    v2_xz.append(ligne[i][8])
    v_1v_xz=[]
    v_2v_xz=[]
    v_o_xz=[]
    for i in range(len(v1_xz)):
        v_1v_xz.append(chiffre(v1_xz[i]))
    for i in range(len(v2_xz)):
        v_2v_xz.append(chiffre(v2_xz[i]))
    for i in range(len(o_xz)):
        v_o_xz.append(chiffre(o_xz[i]))

    bases_xz=[]
    for k in range(len(o_xz)):
        if (verif_xz==1):
            print("---------------------CONFIGURATION XZ-------------------------")
            print("-----------------------Base atome XZ",v_o_xz[k],"-------------")
            print("                                                   ")
            print(ortho(patm[v_o_xz[k]-1],patm[v_1v_xz[k]-1],patm[v_2v_xz[k]-1],13))
            print("                                                   ")
        bases_xz.append([(ortho(patm[v_o_xz[k]-1],patm[v_1v_xz[k]-1],patm[v_2v_xz[k]-1],13))[0],(ortho(patm[v_o_xz[k]-1],patm[v_1v_xz[k]-1],patm[v_2v_xz[k]-1],13))[1],(ortho(patm[v_o_xz[k]-1],patm[v_1v_xz[k]-1],patm[v_2v_xz[k]-1],13))[2]])

    epsilon_t_xz=0

    precision=1e-15
    if (verif_xz==1):
        for t in range(len(bases_xz)):
            epsilon=np.linalg.det(bases_xz[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_xz[t] ," est orthonormee")

##---------------------------------------CONFIGURATION YX -------------------------------#"

    o_yx=[]
    v1_yx=[]
    v2_yx=[]

    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("YX") ):
                o_yx.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_yx.append(residu(ligne[i][4]))
                else:
                    v1_yx.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_yx.append(residu(ligne[i][8]))
                else:
                    v2_yx.append(ligne[i][8])

    v_1v_yx=[]
    v_2v_yx=[]
    v_o_yx=[]
    for i in range(len(v1_yx)):
        v_1v_yx.append(chiffre(v1_yx[i]))
    for i in range(len(v2_yx)):
        v_2v_yx.append(chiffre(v2_yx[i]))
    for i in range(len(o_yx)):
        v_o_yx.append(chiffre(o_yx[i]))



    bases_yx=[]
    for k in range(len(o_yx)):
        if (verif_yx==1):
            print("---------------------CONFIGURATION YX-------------------------")
            print("-----------------------Base atome YX",v_o_yx[k],"-------------")
            print("                                                   ")
            print(ortho(patm[v_o_yx[k]-1],patm[v_1v_yx[k]-1],patm[v_2v_yx[k]-1],21))
            print("                                                   ")
        bases_yx.append([(ortho(patm[v_o_yx[k]-1],patm[v_1v_yx[k]-1],patm[v_2v_yx[k]-1],21))[0],(ortho(patm[v_o_yx[k]-1],patm[v_1v_yx[k]-1],patm[v_2v_yx[k]-1],21))[1],(ortho(patm[v_o_yx[k]-1],patm[v_1v_yx[k]-1],patm[v_2v_yx[k]-1],21
        ))[2]])


    epsilon_t_yx=0

    precision=1e-15
    if (verif_yx==1):
        for t in range(len(bases_yx)):
            epsilon=np.linalg.det(bases_yx[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_yx[t] ," est orthonormee")
    ##------------------------------------------CAS AVEC BISSECTRICE------------------------------------##

    def ortho_biss (atom,voisin1,voisin2):
        u1=[]
        u2=[]
        V1=[]
        V2=[]
        biss=[]
        V1=np.array([voisin1[0]-atom[0],voisin1[1]-atom[1],voisin1[2]-atom[2]])
        V2=np.array([voisin2[0]-atom[0],voisin2[1]-atom[1],voisin2[2]-atom[2]])
        biss=np.array((V1+V2)/2)
        u2=np.array(V2-(np.vdot(V2,biss)/np.vdot(biss,biss))*biss)
        e1=(biss/(np.linalg.norm(biss)))
        e2=((u2/(np.linalg.norm(u2))))
        e3=(np.cross(e1,e2))
        if (A==12):
            return np.array([e1,e2,e3])
        elif (A==13):
            return np.array([e1,e3,e2])
        elif (A==21):
            return np.array([e2,e1,e3])
        elif (A==23):
            return np.array([e3,e1,e2])
        elif (A==31):
            return np.array([e2,e3,e1])
        else:
            return np.array([e3,e2,e1])
    ##--------------------------------------------CONFIGURATION bYX ----------------------------------##
    A=21
    o_byx=[]
    v1_byx=[]
    v2_byx=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("bYX") ):
                o_byx.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_byx.append(residu(ligne[i][4]))
                else:
                    v1_byx.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_byx.append(residu(ligne[i][8]))
                else:
                    v2_byx.append(ligne[i][8])

    v_1v_byx=[]
    v_2v_byx=[]
    v_o_byx=[]
    for i in range(len(v1_byx)):
        v_1v_byx.append(chiffre(v1_byx[i]))
    for i in range(len(v2_yx)):
        v_2v_byx.append(chiffre(v2_byx[i]))
    for i in range(len(o_byx)):
        v_o_byx.append(chiffre(o_byx[i]))


    bases_byx=[]
    for k in range(len(o_byx)):
        if (verif_byx==1):
            print("---------------------CONFIGURATION bYX-------------------------")
            print("-----------------------Base atome bYX",v_o_byx[k],"------------")
            print("                                                   ")
            print((patm[v_o_byx[k]-1],patm[v_1v_byx[k]-1],patm[v_2v_byx[k]-1]))
            print("                                                   ")
        bases_byx.append([(ortho_biss(patm[v_o_byx[k]-1],patm[v_1v_byx[k]-1],patm[v_2v_byx[k]-1]))[0],(ortho_biss(patm[v_o_byx[k]-1],patm[v_1v_byx[k]-1],patm[v_2v_byx[k]-1]))[1],(ortho_biss(patm[v_o_byx[k]-1],patm[v_1v_byx[k]-1],patm[v_2v_byx[k]-1]))[2]])

    epsilon_t_byx=0

    precision=1e-15
    if (verif_byx==1):
        for t in range(len(bases_byx)):
            epsilon=np.linalg.det(bases_byx[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_byx[t] ," est orthonormee")
##--------------------------------------------CONFIGURATION bZX ----------------------------------##
    A=31
    o_bzx=[]
    v1_bzx=[]
    v2_bzx=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("bZX") ):
                o_bzx.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_bzx.append(residu(ligne[i][4]))
                else:
                    v1_bzx.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_bzx.append(residu(ligne[i][8]))
                else:
                    v2_bzx.append(ligne[i][8])

    v_1v_bzx=[]
    v_2v_bzx=[]
    v_o_bzx=[]
    for i in range(len(v1_bzx)):
        v_1v_bzx.append(chiffre(v1_bzx[i]))
    for i in range(len(v2_bzx)):
        v_2v_bzx.append(chiffre(v2_bzx[i]))
    for i in range(len(o_bzx)):
        v_o_bzx.append(chiffre(o_bzx[i]))



    bases_bzx=[]
    for k in range(len(o_bzx)):
        if (verif_bzx==1):
            print("---------------------CONFIGURATION bZX-------------------------")
            print("---------------------Base atome bZX",v_o_bzx[k],"--------------")
            print("                                                   ")
            print(ortho_biss(patm[v_o_bzx[k]-1],patm[v_1v_bzx[k]-1],patm[v_2v_bzx[k]-1]))
            print("                                                   ")
        bases_bzx.append([(ortho_biss(patm[v_o_bzx[k]-1],patm[v_1v_bzx[k]-1],patm[v_2v_bzx[k]-1]))[0],(ortho_biss(patm[v_o_bzx[k]-1],patm[v_1v_bzx[k]-1],patm[v_2v_bzx[k]-1]))[1],(ortho_biss(patm[v_o_bzx[k]-1],patm[v_1v_bzx[k]-1],patm[v_2v_bzx[k]-1]))[2]])


    epsilon_t_bzx=0

    precision=1e-15
    if (verif_bzx==1):
        for t in range(len(bases_bzx)):
            epsilon=np.linalg.det(bases_bzx[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_bzx[t] ," est orthonormee")

##---------------------------------------CONFIGURATION bZY -------------------------------##
    A=32
    o_bzy=[]
    v1_bzy=[]
    v2_bzy=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("bZY") ):
                o_bzy.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_bzy.append(residu(ligne[i][4]))
                else:
                    v1_bzy.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_bzy.append(residu(ligne[i][8]))
                else:
                    v2_bzy.append(ligne[i][8])
    v_1v_bzy=[]
    v_2v_bzy=[]
    v_o_bzy=[]
    for i in range(len(v1_bzy)):
        v_1v_bzy.append(chiffre(v1_bzy[i]))
    for i in range(len(v2_bzy)):
        v_2v_bzy.append(chiffre(v2_bzy[i]))
    for i in range(len(o_bzy)):
        v_o_bzy.append(chiffre(o_bzy[i]))


    bases_bzy=[]
    for k in range(len(o_bzy)):
        if (verif_bzy==1):
            print("-----------------Configuration bZY-----------------")
            print("----------------Base atome",v_o_bzy[k],"-----------")
            print("                                                   ")
            print(ortho_biss(patm[v_o_bzy[k]-1],patm[v_1v_bzy[k]-1],patm[v_2v_bzy[k]-1]))
            print("                                                   ")
        bases_bzy.append([(ortho_biss(patm[v_o_bzy[k]-1],patm[v_1v_bzy[k]-1],patm[v_2v_bzy[k]-1]))[0],(ortho_biss(patm[v_o_bzy[k]-1],patm[v_1v_bzy[k]-1],patm[v_2v_bzy[k]-1]))[1],(ortho_biss(patm[v_o_bzy[k]-1],patm[v_1v_bzy[k]-1],patm[v_2v_bzy[k]-1]))[2]])


    epsilon_t_bzy=0

    precision=0.000000001
    if (verif_bzy==1):
        for t in range(len(bases_bzy)):
            epsilon=np.linalg.det(bases_bzy[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_bzy[t] ," est orthonormee")


##---------------------------------------CONFIGURATION bXY -------------------------------##
    A=12
    o_bxy=[]
    v1_bxy=[]
    v2_bxy=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("bXY") ):
                o_bxy.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_bxy.append(residu(ligne[i][4]))
                else:
                    v1_bxy.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_bxy.append(residu(ligne[i][8]))
                else:
                    v2_bxy.append(ligne[i][8])

    v_1v_bxy=[]
    v_2v_bxy=[]
    v_o_bxy=[]
    for i in range(len(v1_bxy)):
        v_1v_bxy.append(chiffre(v1_bxy[i]))
    for i in range(len(v2_bxy)):
        v_2v_bxy.append(chiffre(v2_bxy[i]))
    for i in range(len(o_bxy)):
        v_o_bxy.append(chiffre(o_bxy[i]))


    bases_bxy=[]
    for k in range(len(o_bxy)):
        if (verif_bxy==1):
            print("------------------------------Configuration bXY-----------")
            print("-----------------------Base atome",v_o_bxy[k],"-----------")
            print("                                                   ")
            print(ortho_biss(patm[v_o_bxy[k]-1],patm[v_1v_bxy[k]-1],patm[v_2v_bxy[k]-1]))
            print("                                                   ")
        bases_bxy.append([(ortho_biss(patm[v_o_bxy[k]-1],patm[v_1v_bxy[k]-1],patm[v_2v_bxy[k]-1]))[0],(ortho_biss(patm[v_o_bxy[k]-1],patm[v_1v_bxy[k]-1],patm[v_2v_bxy[k]-1]))[1],(ortho_biss(patm[v_o_bxy[k]-1],patm[v_1v_bxy[k]-1],patm[v_2v_bxy[k]-1]))[2]])


    epsilon_t_bxy=0

    precision=0.000000001
    if (verif_bxy==1):
        for t in range(len(bases_bxy)):
            epsilon=np.linalg.det(bases_bxy[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_bxy[t] ," est orthonormee")

##---------------------------------------CONFIGURATION bXZ -------------------------------##
    A=13
    o_bxz=[]
    v1_bxz=[]
    v2_bxz=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("bXZ") ):
                o_bxz.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_bxz.append(residu(ligne[i][4]))
                else:
                    v1_bxz.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_bxz.append(residu(ligne[i][8]))
                else:
                    v2_bxz.append(ligne[i][8])

    v_1v_bxz=[]
    v_2v_bxz=[]
    v_o_bxz=[]
    for i in range(len(v1_bxz)):
        v_1v_bxz.append(chiffre(v1_bxz[i]))
    for i in range(len(v2_bxz)):
        v_2v_bxz.append(chiffre(v2_bxz[i]))
    for i in range(len(o_bxz)):
        v_o_bxz.append(chiffre(o_bxz[i]))


    bases_bxz=[]
    for k in range(len(o_bxz)):
        if (verif_bxz==1):
            print("-----------------------Configuration bXZ------------------")
            print("-----------------------Base atome",v_o_bxz[k],"-----------")
            print("                                                   ")
            print(ortho_biss(patm[v_o_bxz[k]-1],patm[v_1v_bxz[k]-1],patm[v_2v_bxz[k]-1]))
            print("                                                   ")
        bases_bxz.append([(ortho_biss(patm[v_o_bxz[k]-1],patm[v_1v_bxz[k]-1],patm[v_2v_bxz[k]-1]))[0],(ortho_biss(patm[v_o_bxz[k]-1],patm[v_1v_bxz[k]-1],patm[v_2v_bxz[k]-1]))[1],(ortho_biss(patm[v_o_bxz[k]-1],patm[v_1v_bxz[k]-1],patm[v_2v_bxz[k]-1]))[2]])


    epsilon_t_bxz=0
    precision=0.000000001
    if (verif_bxz==1):
        for t in range(len(bases_bxz)):
            epsilon=np.linalg.det(bases_bxz[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_bxz[t] ," est orthonormee")

##---------------------------------------CONFIGURATION bYZ -------------------------------##
    A=23
    o_byz=[]
    v1_byz=[]
    v2_byz=[]
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("bYZ") ):
                o_byz.append(ligne[i-1][1])
                if(len(ligne[i][4])>3):
                    v1_byz.append(residu(ligne[i][4]))
                else:
                    v1_byz.append(ligne[i][4])

                if(len(ligne[i][8])>3):
                    v2_byz.append(residu(ligne[i][8]))
                else:
                    v2_byz.append(ligne[i][8])
    v_1v_byz=[]
    v_2v_byz=[]
    v_o_byz=[]
    for i in range(len(v1_byz)):
        v_1v_byz.append(chiffre(v1_byz[i]))
    for i in range(len(v2_byz)):
        v_2v_byz.append(chiffre(v2_byz[i]))
    for i in range(len(o_byz)):
        v_o_byz.append(chiffre(o_byz[i]))

    bases_byz=[]
    for k in range(len(o_byz)):
        if (verif_byz==1):
            print("---------------------Configuration bYZ-----------------")
            print("---------------Base atome",v_o_byz[k],"----------------")
            print("                                                   ")
            print(ortho_biss(patm[v_o_byz[k]-1],patm[v_1v_byz[k]-1],patm[v_2v_byz[k]-1]))
            print("                                                   ")
        bases_byz.append([(ortho_biss(patm[v_o_byz[k]-1],patm[v_1v_byz[k]-1],patm[v_2v_byz[k]-1]))[0],(ortho_biss(patm[v_o_byz[k]-1],patm[v_1v_byz[k]-1],patm[v_2v_byz[k]-1]))[1],(ortho_biss(patm[v_o_byz[k]-1],patm[v_1v_byz[k]-1],patm[v_2v_byz[k]-1]))[2]])



    epsilon_t_byz=0

    precision=0.000000001
    if (verif_byz==1):
        for t in range(len(bases_byz)):
            epsilon=np.linalg.det(bases_byz[t])
            if (abs(epsilon-1.00000000000) <= precision):
                print(epsilon)
                print("la base de latome", v_o_byz[t] ," est orthonormee")



    ##-----------------------------CREATION DES TENSEURS------------------------------------------##
    #variable tempraires
    g=[]
    f=[]
    h=[]
    k=[]
    l=[]
    m=[]
    #--------------------------SECTEUR XY-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("XY")):
                g.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_xy=np.array(g,dtype=float)
    #--------------------------SECTEUR ZX-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("ZX")):
                f.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_zx=np.array(f,dtype=float)
    #--------------------------SECTEUR ZY-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("ZY")):
                h.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_zy=np.array(h,dtype=float)
    #--------------------------SECTEUR YX-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("YX")):
                k.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_yx=np.array(k,dtype=float)
    #--------------------------SECTEUR YZ-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("YZ")):
                l.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_yz=np.array(l,dtype=float)
    #--------------------------SECTEUR XZ-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("XZ")):
                m.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_xz=np.array(m,dtype=float)
    #variables temporaires
    gb=[]
    fb=[]
    hb=[]
    kb=[]
    lb=[]
    mb=[]
    #--------------------------SECTEUR bXY-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("bXY")):
                gb.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_bxy=np.array(gb,dtype=float)
                L=np.reshape(pol_bxy[0],(3,3))
    #--------------------------SECTEUR ZX-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI")  and ligne[i-1][0].startswith("bZX")):
                fb.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_bzx=np.array(fb,dtype=float)
    #--------------------------SECTEUR ZY-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("bZY")):
                hb.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_bzy=np.array(hb,dtype=float)
    #--------------------------SECTEUR YX-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("bYX")):
                kb.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_byx=np.array(kb,dtype=float)
    #--------------------------SECTEUR YZ-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("bYZ")):
                lb.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_byz=np.array(lb,dtype=float)
    #--------------------------SECTEUR XZ-------------------------------------#
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI") and ligne[i-1][0].startswith("bXZ")):
                mb.append([[ligne[i][1],ligne[i][4],ligne[i][5]],[ligne[i][4],ligne[i][2],ligne[i][6]],[ligne[i][5],ligne[i][6],ligne[i][3]]])
                pol_bxz=np.array(mb,dtype=float)

##-----------------------------------MATRICE DE PASSAGE--------------------------------------------#
    P_zx=[]
    for i in range(len(bases_zx)):
        P_zx.append(np.transpose(bases_zx[i]))
    P_xz=[]
    for i in range(len(bases_xz)):
        P_xz.append(np.transpose(bases_xz[i]))
    P_xy=[]
    for i in range(len(bases_xy)):
        P_xy.append(np.transpose(bases_xy[i]))
    P_yx=[]
    for i in range(len(bases_yx)):
        P_yx.append(np.transpose(bases_yx[i]))
    P_zy=[]
    for i in range(len(bases_zy)):
        P_zy.append(np.transpose(bases_zy[i]))
    P_yz=[]
    for i in range(len(bases_yz)):
        P_yz.append(np.transpose(bases_yz[i]))
#------------------------cas avec bissectrice -----------------------------------------------------#
    P_bzx=[]
    for i in range(len(bases_bzx)):
        P_bzx.append(np.transpose(bases_bzx[i]))
    P_bxz=[]
    for i in range(len(bases_bxz)):
        P_bxz.append(np.transpose(bases_bxz[i]))
    P_bxy=[]
    for i in range(len(bases_bxy)):
        P_bxy.append(np.transpose(bases_bxy[i]))
    P_byx=[]
    for i in range(len(bases_byx)):
        P_byx.append(np.transpose(bases_byx[i]))
    P_bzy=[]
    for i in range(len(bases_bzy)):
        P_bzy.append(np.transpose(bases_bzy[i]))
    P_byz=[]
    for i in range(len(bases_byz)):
        P_byz.append(np.transpose(bases_byz[i]))


##------- CHANGEMENT DE BASE DES TENSEURS DE B(molecule) vers B(atom) -----------##
    #-------------------------Cas XY------------------------------------------
    pol_t_xy=[]

    for i in range(len(P_xy)):
        pol_t_xy.append(np.dot(np.linalg.inv(P_xy[i]),np.dot(pol_xy[i],P_xy[i])))
        if (verif_xy==1):
            print("----------------------------Configuration XY--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_xy[i],"--------------")
            print(pol_t_xy[i])
            print("                                                                           ")
            print("vp de pol_zx[0] avant changement",np.linalg.eig((pol_yx[i]).reshape(3,3))[0])
            print("vp de pol_t_zx ",np.linalg.eig((pol_t_xy[i]).reshape(3,3))[0])

    #-------------------------Cas XZ-----------------------------------------------
    pol_t_xz=[]
    for i in range(len(P_xz)):
        pol_t_xz.append(np.dot(np.linalg.inv(P_xz[i]),np.dot(pol_xz[i],P_xz[i])))
        if (verif_xz==1):
            print("----------------------------Configuration XZ--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_xz[i],"--------------")
            print(pol_t_xz[i])
            print("                                                                           ")
            print("vp de pol_xz[0] avant changement",np.linalg.eig((pol_xz[i]).reshape(3,3))[0])
            print("vp de pol_t_xz ",np.linalg.eig((pol_t_xz[i]).reshape(3,3))[0])
    #-------------------------Cas YX-----------------------------------------------------
    pol_t_yx=[]
    for i in range(len(P_yx)):
        pol_t_yx.append(np.dot(np.linalg.inv(P_yx[i]),np.dot(pol_yx[i],P_yx[i])))
        if (verif_yx==1):
            print("----------------------------Configuration YX--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_yx[i],"--------------")
            print(pol_t_yx[i])
            print("                                                                           ")
            print("vp de pol_yx[0] avant changement",np.linalg.eig((pol_yx[i]).reshape(3,3))[0])
            print("vp de pol_t_yx ",np.linalg.eig((pol_t_yx[i]).reshape(3,3))[0])
    #-------------------------Cas YZ----------------------------------------------------
    pol_t_yz=[]
    for i in range(len(P_yz)):
        pol_t_yz.append(np.dot(np.linalg.inv(P_yz[i]),np.dot(pol_yz[i],P_yz[i])))
        if (verif_yz==1):
            print("----------------------------Configuration YZ--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_yz[i],"--------------")
            print(pol_t_yz[i])
            print("                                                                           ")
            print("vp de pol_yz[0] avant changement",np.linalg.eig((pol_yz[i]).reshape(3,3))[0])
            print("vp de pol_t_yz ",np.linalg.eig((pol_t_xz[i]).reshape(3,3))[0])
    #--------------------------Cas ZX-----------------------------------------------
    pol_t_zx=[]
    for i in range(len(P_zx)):
        pol_t_zx.append(np.dot(np.linalg.inv(P_zx[i]),np.dot(pol_zx[i],P_zx[i])))
        if (verif_zx==1):
            print("----------------------------Configuration ZX--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_zx[i],"--------------")
            print(pol_t_zx[i])
            print("                                                                           ")
            print("vp de pol_zx avant changement",np.linalg.eig((pol_zx[i]).reshape(3,3))[0])
            print("vp de pol_t_zx ",np.linalg.eig((pol_t_zx[i]).reshape(3,3))[0])

    #--------------------------Cas ZY-----------------------------------------------
    pol_t_zy=[]
    for i in range(len(P_zy)):
        pol_t_zy.append(np.dot(np.linalg.inv(P_zy[i]),np.dot(pol_zy[i],P_zy[i])))
        if (verif_zy==1):
            print("----------------------------Configuration ZY--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_zy[i],"--------------")
            print(pol_t_zy[i])
            print("                                                                           ")
            print("vp de pol_zy[0] avant changement",np.linalg.eig((pol_zy[i]).reshape(3,3))[0])
            print("vp de pol_t_zy ",np.linalg.eig((pol_t_zy[i]).reshape(3,3))[0])

    #-------------------------Cas bXY------------------------------------------#
    pol_t_bxy=[]
    for i in range(len(P_bxy)):
        pol_t_bxy.append(np.dot(np.linalg.inv(P_bxy[i]),np.dot(pol_bxy[i],P_bxy[i])))
        if (verif_bxy==1):
            print("----------------------------Configuration bxy--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_bxy[i],"--------------")
            print(pol_t_bxy[i])
            print("                                                                           ")
            print("vp de pol_bxy[0] avant changement",np.linalg.eig((pol_bxy[i]).reshape(3,3))[0])
            print("vp de pol_t_bxy ",np.linalg.eig((pol_t_bxy[i]).reshape(3,3))[0])
    #-------------------------Cas bXZ-----------------------------------------------
    pol_t_bxz=[]
    for i in range(len(P_bxz)):
        pol_t_bxz.append(np.dot(np.linalg.inv(P_bxz[i]),np.dot(pol_bxz[i],P_bxz[i])))
        if (verif_bxz==1):
            print("----------------------------Configuration bxz--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_bxz[i],"--------------")
            print(pol_t_bxz[i])
            print("                                                                           ")
            print("vp de pol_bxz[0] avant changement",np.linalg.eig((pol_bxz[i]).reshape(3,3))[0])
            print("vp de pol_t_bxz ",np.linalg.eig((pol_t_bxz[i]).reshape(3,3))[0])
    #-------------------------Cas bYX-----------------------------------------------------
    pol_t_byx=[]
    for i in range(len(P_byx)):
        pol_t_byx.append(np.dot(np.linalg.inv(P_byx[i]),np.dot(pol_byx[i],P_byx[i])))
        if (verif_byx==1):
            print("----------------------------Configuration byx--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_byx[i],"--------------")
            print(pol_t_byx[i])
            print("                                                                           ")
            print("vp de pol_byx[0] avant changement",np.linalg.eig((pol_byx[i]).reshape(3,3))[0])
            print("vp de pol_t_byx ",np.linalg.eig((pol_t_byx[i]).reshape(3,3))[0])
    #-------------------------Cas bYZ----------------------------------------------------
    pol_t_byz=[]
    for i in range(len(P_byz)):
        pol_t_byz.append(np.dot(np.linalg.inv(P_byz[i]),np.dot(pol_byz[i],P_byz[i])))
        if (verif_byz==1):
            print("----------------------------Configuration byz--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_byz[i],"--------------")
            print(pol_t_byz[i])
            print("                                                                           ")
            print("vp de pol_byz[0] avant changement",np.linalg.eig((pol_yz[i]).reshape(3,3))[0])
            print("vp de pol_t_byz ",np.linalg.eig((pol_t_byz[i]).reshape(3,3))[0])
    #--------------------------Cas bZX-----------------------------------------------
    pol_t_bzx=[]
    for i in range(len(P_bzx)):
        pol_t_bzx.append(np.dot(np.linalg.inv(P_bzx[i]),np.dot(pol_bzx[i],P_bzx[i])))
        if (verif_bzx==1):
            print("----------------------------Configuration bzx--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_bzx[i],"--------------")
            print(pol_t_bzx[i])
            print("                                                                           ")
            print("vp de pol_bzx[0] avant changement",np.linalg.eig((pol_bzx[i]).reshape(3,3))[0])
            print("vp de pol_t_bzx ",np.linalg.eig((pol_t_bzx[i]).reshape(3,3))[0])
    #--------------------------Cas bZY-----------------------------------------------
    pol_t_bzy=[]
    for i in range(len(P_bzy)):
        pol_t_bzy.append(np.dot(np.linalg.inv(P_bzy[i]),np.dot(pol_bzy[i],P_bzy[i])))
        if (verif_bzy==1):
            print("----------------------------Configuration bzy--------------------------")
            print("Le tenseur de polarisabilité associé à l'atome",v_o_bzy[i],"--------------")
            print(pol_t_bzy[i])
            print("                                                                           ")
            print("vp de pol_bzy[0] avant changement",np.linalg.eig((pol_bzy[i]).reshape(3,3))[0])
            print("vp de pol_t_bzy ",np.linalg.eig((pol_t_bzy[i]).reshape(3,3))[0])



    ##ecriture des données dans un txt

    type_atom=[]
    num_atom=[]
    name_atom=[]
    val_pol=[]
    config_atom=[]
    residu_atom=[]
    voisins_atom=[]
    #recuperation des coefficients des dipoles et quadripoles
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j].startswith("UANI")):
                if ligne [i-1][12]=="DIP":
                    val_pol.append([ligne[i+1][2],ligne[i+1][3],ligne[i+1][4],'0','0','0','0','0'])
                else:
                    val_pol.append([ligne[i+1][2],ligne[i+1][3],ligne[i+1][4],ligne[i+1][5],ligne[i+1][6],ligne[i+1][7],ligne[i+1][8],ligne[i+1][9]])


    #recupération des infos utiles a reinscrire en OUTPUT
    for i in range(len(ligne)):
        for j in range(len(ligne[i])):
            if (ligne[i][j]=="ATOM"):
                type_atom.append(ligne[i][18])
                num_atom.append(int(ligne[i][1]))
                name_atom.append(ligne[i][2])
                config_atom.append(ligne[i+1][0])
                residu_atom.append(ligne[i][4])
                voisins_atom.append([ligne[i+1][4],ligne[i+1][8]])


    #trouver le type atomique d'un atome en fontion de son numéro
    def find_type_atom(numero):
        return type_atom[numero-1]
    #trouver le nom d'un atome en fonctione de son numéro
    def find_name_atom(numero):
        return name_atom[int(numero)-1]
    #variables regroupants tout les atomes en fonction de leur config
    all_V=[]
    #contient les listes des polarizations de chaque config
    all_pol=[]
    all_V=[v_o_xy,v_o_yx,v_o_yz,v_o_zy,v_o_zx,v_o_xz,v_o_bxy,v_o_byx,v_o_byz,v_o_bzy,v_o_bzx,v_o_bxz]
    all_pol=[pol_t_xy,pol_t_yx,pol_t_yz,pol_t_zy,pol_t_zx,pol_t_xz,pol_t_bxy,pol_t_byx,pol_t_byz,pol_t_bzy,pol_t_bzx,pol_t_bxz]
    #contient toutes les listes de bases de chaque config
    all_bases=[bases_xy,bases_yx,bases_yz,bases_zy,bases_zx,bases_bxz,bases_bxy,bases_byx,bases_byz,bases_bzy,bases_bzx,bases_bxz]
    def find_pol(numero):
            for i in range(12):
                for j in range(len(all_V[i])):
                    if (all_V[i][j]== numero) :
                        return all_pol[i][j]
    def find_bases(numero):
            for i in range(12):
                for j in range(len(all_V[i])):
                    if (all_V[i][j]== numero) :
                        return all_bases[i][j]


    ##----------------------Sorti des données de tout les atomes du set en .txt-------------------------------------------##

    os.chdir(path_to_txtfiles)       #repertoire darrivee du fichier OUTPUT.txt en sortie
    extfile=path_to_txtfiles
    if (x==0):
        with open ("pol_OUTPUT_TRANSFORM.txt", "w") as f:
            f.write(" format : N° NAMEATOM MOL RESIDU CONFIG VOIS1 VOIS2 XX YY ZZ XY XZ YZ DIP1 DIP2 DIP3 QUA1 QUA2 QUA3 QUA4 QUA5 NAMEFILE"'\n')
            f.write('\n')
            f.write('\n')

        with open ("pol_OUTPUT_TRANSFORM.txt", "a") as f:

            for i in range(len(num_atom)):
                f.write ( str(num_atom[i])+'\t'+ (find_name_atom(num_atom[i]))+'\t'+"MOL"+'\t'+residu_atom[i]+'\t'+config_atom[i]+'\t'+voisins_atom[i][0]+'\t'+voisins_atom[i][1]+'\t'+str(find_pol(num_atom[i])[0][0])+'\t'+str(find_pol(num_atom[i])[1][1])+'\t'+ str(find_pol(num_atom[i])[2][2])+'\t'+ str(find_pol((num_atom[i]))[1][0])+'\t'+ str(find_pol((num_atom[i]))[2][0]) + '\t'+ str(find_pol(num_atom[i])[2][1]) + '\t'+find_type_atom(num_atom[i])+'\t'+val_pol[i][0]+'\t'+val_pol[i][1]+'\t'+val_pol[i][2]+'\t'+val_pol[i][3]+'\t'+val_pol[i][4]+'\t'+val_pol[i][5]+'\t'+val_pol[i][6]+'\t'+val_pol[i][7]+'\t'+str(find_bases(num_atom[i])[0][0])+'\t'+str(find_bases(num_atom[i])[0][1])+'\t'+str(find_bases(num_atom[i])[0][2])+'\t'+str(find_bases(num_atom[i])[1][0])+'\t'+str(find_bases(num_atom[i])[1][1])+'\t'+str(find_bases(num_atom[i])[1][2])+'\t'+str(find_bases(num_atom[i])[2][0])+'\t'+str(find_bases(num_atom[i])[2][1])+'\t'+str(find_bases(num_atom[i])[2][2])+'\t'+fichiers[x] +'\n')
    else:
        with open ("pol_OUTPUT_TRANSFORM.txt", "a") as f:

            for i in range(len(num_atom)):
                f.write ( str(num_atom[i])+'\t'+ (find_name_atom(num_atom[i]))+'\t'+"MOL"+'\t'+residu_atom[i]+'\t'+config_atom[i]+'\t'+voisins_atom[i][0]+'\t'+voisins_atom[i][1]+'\t'+str(find_pol(num_atom[i])[0][0])+'\t'+str(find_pol(num_atom[i])[1][1])+'\t'+ str(find_pol(num_atom[i])[2][2])+'\t'+ str(find_pol((num_atom[i]))[1][0])+'\t'+ str(find_pol((num_atom[i]))[2][0]) + '\t'+ str(find_pol(num_atom[i])[2][1]) + '\t'+find_type_atom(num_atom[i])+'\t'+val_pol[i][0]+'\t'+val_pol[i][1]+'\t'+val_pol[i][2]+'\t'+val_pol[i][3]+'\t'+val_pol[i][4]+'\t'+val_pol[i][5]+'\t'+val_pol[i][6]+'\t'+val_pol[i][7]+'\t'+str(find_bases(num_atom[i])[0][0])+'\t'+str(find_bases(num_atom[i])[0][1])+'\t'+str(find_bases(num_atom[i])[0][2])+'\t'+str(find_bases(num_atom[i])[1][0])+'\t'+str(find_bases(num_atom[i])[1][1])+'\t'+str(find_bases(num_atom[i])[1][2])+'\t'+str(find_bases(num_atom[i])[2][0])+'\t'+str(find_bases(num_atom[i])[2][1])+'\t'+str(find_bases(num_atom[i])[2][2])+'\t'+fichiers[x] +'\n')

##-------------------------COPIE DES FICHIERS ORIGINAUX DU BIG SET EN AJOUTANT LES "TPOL"----------------------------------#
    ligne3=[]
    ligne4=[]
    os.chdir(path_to_the_bigset)
    with open(fichiers[x], 'r') as file:
        data=file.readlines()
        data2=data
        for li in data:  #fichier texte dans une liste
            s = li.strip ("\n\r")
            l = s .split (" ")
            ligne4.append (l) #on découpe en colonne pour pouvoir choisir les mots
            ligne3.append (l)


    q=0 #incrémentation
    if tpol=='ON':
#On fait une première boucle sur ligne4 car on ne peut pas ajouter de ligne sur une liste que l'on itère
        for i in range(len(ligne4)):
            for j in range(len(ligne4[i])):
                if ligne4[i][j]=="ATOM":
                    if ligne4[i+1][12]=='OCT':

                        ligne3[i+5]=('TPOL'+'\t'+str(find_pol(int(ligne4[i][1]))[0][0])+'\t'+str(find_pol(int(ligne4[i][1]))[1][1])+'\t'+ str(find_pol(int(ligne4[i][1]))[2][2])+'\t'+str(find_pol((int(ligne4[i][1])))[1][0])+'\t'+str(find_pol((int(ligne4[i][1])))[2][0])+'\t'+str(find_pol(int(ligne4[i][1]))[2][1])+'\n')

                        data2[i+5]=('TPOL'+'\t'+str(find_pol(int(ligne4[i][1]))[0][0])+'\t'+str(find_pol(int(ligne4[i][1]))[1][1])+'\t'+ str(find_pol(int(ligne4[i][1]))[2][2])+'\t'+str(find_pol((int(ligne4[i][1])))[1][0])+'\t'+str(find_pol((int(ligne4[i][1])))[2][0])+'\t'+str(find_pol(int(ligne4[i][1]))[2][1])+'\n')


                        # print('TPOL',q)
                        q=q+1

                    elif ligne4[i+1][12]=='HEX':

                        ligne3[i+6]=('TPOL'+'\t'+str(find_pol(int(ligne4[i][1]))[0][0])+'\t'+str(find_pol(int(ligne4[i][1]))[1][1])+'\t'+ str(find_pol(int(ligne4[i][1]))[2][2])+'\t'+str(find_pol((int(ligne4[i][1])))[1][0])+'\t'+str(find_pol((int(ligne4[i][1])))[2][0])+'\t'+str(find_pol(int(ligne4[i][1]))[2][1])+'\n')

                        data2[i+6]=('TPOL'+'\t'+str(find_pol(int(ligne4[i][1]))[0][0])+'\t'+str(find_pol(int(ligne4[i][1]))[1][1])+'\t'+ str(find_pol(int(ligne4[i][1]))[2][2])+'\t'+str(find_pol((int(ligne4[i][1])))[1][0])+'\t'+str(find_pol((int(ligne4[i][1])))[2][0])+'\t'+str(find_pol(int(ligne4[i][1]))[2][1])+'\n')


                        # print('TPOL',q)
                        q=q+1

                    else:

                        ligne3[i+4]=('TPOL'+'\t'+str(find_pol(int(ligne4[i][1]))[0][0])+'\t'+str(find_pol(int(ligne4[i][1]))[1][1])+'\t'+ str(find_pol(int(ligne4[i][1]))[2][2])+'\t'+str(find_pol((int(ligne4[i][1])))[1][0])+'\t'+str(find_pol((int(ligne4[i][1])))[2][0])+'\t'+str(find_pol(int(ligne4[i][1]))[2][1])+'\n')


                        data2[i+4]=('TPOL'+'\t'+str(find_pol(int(ligne4[i][1]))[0][0])+'\t'+str(find_pol(int(ligne4[i][1]))[1][1])+'\t'+ str(find_pol(int(ligne4[i][1]))[2][2])+'\t'+str(find_pol((int(ligne4[i][1])))[1][0])+'\t'+str(find_pol((int(ligne4[i][1])))[2][0])+'\t'+str(find_pol(int(ligne4[i][1]))[2][1])+'\n')

                        #
                        # print('TPOL',q)
                        q=q+1

    #On compte le nombre de fois qu'apparait TPOL, c'est le nombre de saut de ligne à rajouter
        ft=0
        for i in range(len(data2)):
            if data2[i].startswith("TPOL"):
                ft=ft+1

    #On fait le saut de ligne ( donc ajout ) seulement maintenant
        for i in range(len(data2)+ft):
            if data2[i].startswith("TPOL"):
                data2.insert(i+1,'\n')



        with open('TPOL_'+fichiers[x], 'a') as new:
            new.writelines(data2)
        shutil.move('TPOL_'+fichiers[x], desti)


print(' ')
print('OUTPUT file ready for further work created in :',path_to_txtfiles)
print("Data treated in : --- %s seconds ---" % (time.time() - start_time))
if code_exit=='ON':
    ##---------------------2/Detrmination des types atomiques et obtention de leur différentes polarisabilités -------------------------------##

    ##-----------------------------determination des types atomiques (premier tri avec grande précision)--------------------------------##



    #cette étape nous a donnés un nombre de types atomiques différents = 201, or cela est beaucoup trop précis comme tri, on a raffiné manuellement par comparaison des densités electroniques que l'on pouvait réduire ce nombre à 60
    def lettre(X):
        return str(re.sub('[0123456789_]','',X))

    #repertoire ou vous voulez creer le fichier des types atomiques
    os.chdir(path_to_txtfiles)

    ligne2=[]
    with open ("pol_OUTPUT_TRANSFORM.txt", "r") as file:
            for li in file :
                s = li.strip ("\n""\r""\t")
                l = s.split ('\t')
                ligne2.append (l)

    def str_compare(j,k):
        if j==k:
            return True
        else:
            return False


    #liste de tout les elements du fichier pour chaque ligne
    element=[]
    for i in range(3,len(ligne2)):
        for j in range(len(ligne2[i])):
            if ligne2[i][j]=="MOL":
                element.append(ligne2[i][13])
    #tri et trouve tout les différents type delement chimiques et les assignent dans une liste
    file_element=[]
    for i in range(len(element)):
        if element[i] not in file_element:
            file_element.append(element[i])
        else:
            pass

    #rassemble dans une liste de vecteurs les 8 coefficients dipolaire et quadruolaires de chaque ligne
    pval=[]
    for i in range(3,len(ligne2)):
        pval.append([ligne2[i][14],ligne2[i][15],ligne2[i][16],ligne2[i][17],ligne2[i][18],ligne2[i][19],ligne2[i][20],ligne2[i][21]])

    #toute les configurations possible des 8 coefficients pval dans le .txt
    file_pval=[]
    for i in range(len(pval)):
        if pval[i] not in file_pval:
            file_pval.append(pval[i])
        else:
            pass

    #toute les configurations géométriques possibles dans le fichier
    file_config=['XY','YX','YZ','ZY','ZX','XZ','bXY','bYX','bYZ','bZY','bZX','bXZ']

    #liste des config par ligne
    config=[]
    for i in range(3,len(ligne2)):
            config.append(ligne2[i][4])
    #renvoie lelement chimique du nom du voisin (ex voisin('C23') return C
    def voisin(x):
            return str(re.sub('[0123456789_]' ,'',x))

    #creation des liste de voisin 1 et 2 pour chaque ligne
    voisin1=[]
    voisin2=[]

    for i in range(3,len(ligne2)):
            voisin1.append(voisin(ligne2[i][5]))
            voisin2.append(voisin(ligne2[i][6]))
    #nom du fichier
    filename=[]
    for i in range(3,len(ligne2)):
        filename.append(ligne2[i][31])

    #nom atome
    name=[]
    for i in range(3,len(ligne2)):
        name.append(ligne2[i][0])

    #construction du vecteur [pval , element , config , v1 , v2 ] pour chaque ligne
    type_atom=[]
    for i in range(len(pval)):
            type_atom.append([pval[i],element[i],config[i],voisin1[i],voisin2[i]])


    #recherche des types atomiques dans le fichier
    a=0
    file_type_atom=[]
    for i in range(len(pval)):
            if ([pval[i],element[i],config[i],voisin1[i],voisin2[i]] not in file_type_atom):
                    file_type_atom.append([pval[i],element[i],config[i],voisin1[i],voisin2[i]])
                    a=a+1


    #verification si les types atomiques ont bien uniques(ici les types sont brut, on va par la suite extraires ces types et les analyser)
    for i in range(len(file_type_atom)):
            for j in range(len(file_type_atom)):
                    if (str_compare(file_type_atom[i],file_type_atom[j])=='True'and i!=j):
                            print('ca marche pas',i,j)
                    else:
                            pass


    type_atom_infos=[]
    nb_occur=[]
    #vecteur contenant le nombre d'occurence dun type atomique et ses parametres
    for i in range(len(file_type_atom)):
            type_atom_infos.append([type_atom.count(file_type_atom[i]),file_type_atom[i]])
            nb_occur.append(type_atom.count(file_type_atom[i]))

    type_atom_inf=[]
    type_atom_unique=[]
    type_atom_sup=[]
    for j in range(len(type_atom_infos)):
            if (type_atom_infos[j][0] >=2 and type_atom_infos[j][0] <= 10):
                    type_atom_inf.append(type_atom_infos[j][1])
            elif (type_atom_infos[j][0] == 1):
                    type_atom_unique.append(type_atom_infos[j][1])
            else:
                    type_atom_sup.append(type_atom_infos[j][1])

    #création du fichier du premier tri grossier (201 types)
    os.chdir(path_to_txtfiles)
    with open("Type_atom_output.txt",'w') as f:
        f.write('Nombre de type atomiques différents dans le BigSet:'+'\t'+str(len(file_type_atom))+'\n')
        f.write('Nombre de type atomiques à occurence unique dans le BigSet:'+'\t'+str(len(type_atom_unique))+'\n')
        f.write('Nombre de type atomiques à plus de 10 occurences dans le BigSet:'+'\t'+str(len(type_atom_sup))+'\n')
        f.write('Nombre de type atomiques avec un nombre doccurence compris entre 2 et 10 dans le BigSet:'+'\t'+str(len(type_atom_inf))+'\n')
        f.write('\n')
        f.write('Liste des types atomiques'+'\n')
        f.write('\n')
        for i in range(len(type_atom_infos)):
                    f.write(str(type_atom_infos[i][0])+'\t'+type_atom_infos[i][1][0][0]+'\t'+type_atom_infos[i][1][0][1]+'\t'+type_atom_infos[i][1][0][2]+'\t'+type_atom_infos[i][1][0][3]+'\t'+type_atom_infos[i][1][0][4]+'\t'+type_atom_infos[i][1][0][5]+'\t'+type_atom_infos[i][1][0][6]+'\t'+type_atom_infos[i][1][0][7]+'\t')
                    f.write(type_atom_infos[i][1][1]+'\t'+type_atom_infos[i][1][2]+'\t'+type_atom_infos[i][1][3]+'\t'+type_atom_infos[i][1][4]+'\n')


    ##affinage du nombre de types atomiques ##

    ##Fonction utile pour le tri atomique (multipole to fichiers/atom)
    if (find_type_atom_by_dip_qua_config =='ON'):
        print('entrer un type atomique')
        type_precise= [str(x) for x in input().split()]
        for i in range(len(type_atom)):
            if(type_precise==[type_atom[i][0][0],type_atom[i][0][1],type_atom[i][0][2],type_atom[i][0][3],type_atom[i][0][4],type_atom[i][0][5],type_atom[i][0][6],type_atom[i][0][7],type_atom[i][1],type_atom[i][2],type_atom[i][3],type_atom[i][4]]):
                print(filename[i],print(name[i]))

    #sorti du fichier txt des 201 types atomiques "brut" pour affinage du nombre de types sur excel par la suite a la main.
    new_tri=[]
    for i in range(len(file_type_atom)):
        if ([file_type_atom[i][1],file_type_atom[i][2],file_type_atom[i][3],file_type_atom[i][4]] not in new_tri):
            new_tri.append([file_type_atom[i][1],file_type_atom[i][2],file_type_atom[i][3],file_type_atom[i][4]])

    new_tri_dq=[]
    for i in range(len(file_type_atom)):
        for j in range(len(new_tri)):
            if ([file_type_atom[i][1],file_type_atom[i][2],file_type_atom[i][3],file_type_atom[i][4]]==new_tri[j]):
                new_tri_dq.append([file_type_atom[i][0],j])

    aa=sorted(new_tri_dq, key=itemgetter(1))


    os.chdir(path_to_txtfiles)
    with open("TRI_TYPE_ATOM.txt",'w') as f:
        f.write('Tri par association config v1 v2 et element'+'\n')
        f.write('pour chaque : O bxy H H , on a une liste de differentes association de dipoles quadrupoles'+'\n')
        f.write('\n')
        for i in range(len(new_tri)):
            f.write(new_tri[i][0]+'\t'+new_tri[i][1]+'\t'+new_tri[i][2]+'\t'+new_tri[i][3]+'\t'+str(i)+'\n')
        f.write('\n')
        for i in range(len(aa)):
            f.write(aa[i][0][0]+'\t'+aa[i][0][1]+'\t'+aa[i][0][2]+'\t'+aa[i][0][3]+'\t'+aa[i][0][4]+'\t'+aa[i][0][5]+'\t'+aa[i][0][6]+'\t'+aa[i][0][7]+'\t'+str(aa[i][1])+'\n')


    ##-------------------------- TRI DES POLARISATIONS ATOMIQUES PAR TYPE ATOMIQUE ------------------------##



    os.chdir(path_to_txtfiles)

    ligne5=[]
    with open (path_to_liste, "r") as file:
            for li in file :
                s = li.strip ("\n""\r""\t")
                l = s.split ('\t')
                ligne5.append (l)

    #declaration des variables prealables, celles ci sont propres au set utilisé
    #rechiffrage du code .par00 en lettre (utilisé pour laffinage du tri des types atomiques sur excel ),(ex: num[3]=3 --> codepar[3]=['H','ZX','O','C'] : le code pour la configuration ['H','ZX','O','C'] est le numéro 3
    num=[]
    for i in range(40):
        num.append(i)

    numéro=list(range(0,40))
    numéro_type=[]
    numéro_type.append([['O','bXY','H','H'],
    ['H','ZX','O','H'],
    ['O','bXY','C','H'],
    ['H','ZX','O','C'],
    ['C','bZX','O','H'],
    ['H','ZX','C','O'],
    ['N','ZX','C','H'],
    ['H','ZX','N','C'],
    ['C','ZX','N','H'],
    ['H','ZX','C','N'],
    ['C','ZX','C','H'],
    ['H','ZX','C','C'],
    ['C','bXY','O','N'],
    ['O','XY','C','N'],
    ['N','bXY','C','C'],
    ['C','bXY','C','C'],
    ['C','bXY','N','C'],
    ['C','bXY','N','N'],
    ['C','bXY','O','O'],
    ['O','XY','C','O'],
    ['N','XY','C','H'],
    ['C','bXY','H','H'],
    ['C','ZX','C','C'],
    ['C','bZX','C','C'],
    ['H','ZX','C','H'],
    ['C','ZX','C','O'],
    ['O','bXY','C','C'],
    ['C','bZX','O','C'],
    ['Cl','bZX','C','N'],
    ['C','bZX','S','C'],
    ['H','ZX','C','S'],
    ['S','bXY','C','C'],
    ['C','ZX','S','H'],
    ['C','bZX','N','C'],
    ['N','ZX','H','H'],
    ['H','ZX','N','H'],
    ['H','ZX','O','P'],
    ['O','XY','P','O'],
    ['O','bXY','P','H'],
    ['P','bZX','O','O']])



    import copy
    ligne6 = copy.deepcopy(ligne5)

    for i in range(len(ligne6)):
        for j in range(len(numéro)):
            if int(ligne6[i][8])==numéro[j]:

                ligne5[i].remove(ligne5[i][8])
                ligne5[i].append(numéro_type[0][j][0])
                ligne5[i].append(numéro_type[0][j][1])
                ligne5[i].append(numéro_type[0][j][2])
                ligne5[i].append(numéro_type[0][j][3])
    ##--------------------------------------------------------------------------------------------------------------#
    #chagement de variable car type_atom a été construit d'une manière peu pratique précédement, ALL_ATOM contient les dipoles, quadripoles ainsi que la configuration, le nom de l'atom et les 2 voisins
    ALL_ATOM=[]
    for i in range(len(type_atom)):
        ALL_ATOM.append([type_atom[i][0][0],type_atom[i][0][1],type_atom[i][0][2],type_atom[i][0][3],type_atom[i][0][4],type_atom[i][0][5],type_atom[i][0][6],type_atom[i][0][7],type_atom[i][1],type_atom[i][2],type_atom[i][3],type_atom[i][4]])


    #création dun dictionnaire pour les noms des 60 types atomiques présents dans le ficier "liste_type_atom_done.txt" qui nous permetrrons de trier ceux du fichier ou il y a tout les atomes
    name_type_final=[]
    for i in range(len(ligne5)):
        if (ligne5[i][8] not in name_type_final):
            name_type_final.append(ligne5[i][8])
    d = {}
    for i in range(len(name_type_final)):
        d["type"+name_type_final[i]] = []


    for j in range(len(name_type_final)):
        for i in range(len(ligne5)):
            if (ligne5[i][0]!='' and ligne5[i][8]==name_type_final[j]):
                d["type"+name_type_final[j]].append([ligne5[i][0],ligne5[i][1],ligne5[i][2],ligne5[i][3],ligne5[i][4],ligne5[i][5],ligne5[i][6],ligne5[i][7],ligne5[i][9],ligne5[i][10],ligne5[i][11],ligne5[i][12]])

    #dictionnaire des polarisations
    dpol= {}
    for i in range(len(name_type_final)):
        dpol["pol" + name_type_final[i] ] = []


    #remplissage du dictionnaire des polarisations, pour afficher toutes les polarisations d'un type : dpol["type"+nomdutype] ex : dpol["typeC209"].

    for i in range(3,len(ligne2)):
        for j in range(len(d)):
            if (ALL_ATOM[i-3] in d["type"+name_type_final[j]]):
                dpol["pol" + name_type_final[j]].append([ligne2[i][7],ligne2[i][8],ligne2[i][9],ligne2[i][10],ligne2[i][11],ligne2[i][12]])



    ##-------------PROCEDES STATISTIQUES DANALYSE DES POL-------------------------##

    #MOYENNES DES VALEURS PROPRES DES TPOL (travail dans le fichier txt "POL_TYPE_ATOM")

    os.chdir(path_to_txtfiles)
    #dictionnaire des valeurs propres et vecteurs propre pour chaque tenseur
    eigv_pol= {}
    for i in range(len(name_type_final)):
        eigv_pol["eigv_pol" + name_type_final[i] ] = []


    for i in range(len(dpol)):
        for j in range(len(dpol["pol"+name_type_final[i]])):
            t=np.array([[dpol["pol"+name_type_final[i]][j][0],dpol["pol"+name_type_final[i]][j][3],dpol["pol"+name_type_final[i]][j][4]],[dpol["pol"+name_type_final[i]][j][3],dpol["pol"+name_type_final[i]][j][1],dpol["pol"+name_type_final[i]][j][5]],[dpol["pol"+name_type_final[i]][j][4],dpol["pol"+name_type_final[i]][j][5],dpol["pol"+name_type_final[i]][j][2]]],dtype=float)
            b=np.linalg.eigh(t)
            eigv_pol["eigv_pol"+name_type_final[i]].append(b)

    #dictionnaire des 3 vp
    alleig={}
    for i in range(len(name_type_final)):
        alleig[name_type_final[i]]=[]

    #dictionnaire des ecarts types
    std={}
    for i in range(len(name_type_final)):
        std[name_type_final[i]]=[]

    #tracé des histogramme de dispersion pour chaque type
    for i in range(len(eigv_pol)):
        vx=[]
        vy=[]
        vz=[]
        for j in range(len(eigv_pol["eigv_pol"+name_type_final[i]])):
            vx.append(eigv_pol["eigv_pol"+name_type_final[i]][j][0][0])
            vy.append(eigv_pol["eigv_pol"+name_type_final[i]][j][0][1])
            vz.append(eigv_pol["eigv_pol"+name_type_final[i]][j][0][2])
            Vx=np.round(vx,7)
            Vy=np.round(vy,7)
            Vz=np.round(vz,7)
            Vx.sort()
            Vy.sort()
            Vz.sort()
        alleig[name_type_final[i]].append([Vx,Vy,Vz])

    def find_distr(type):
        sns.set_style("whitegrid")
        findhist=str(type)
        g=sns.displot(data=alleig[findhist][0][0],bins=len(str(len(alleig[name_type_final[q]][0][0])))*20,kde=True,color='red')
        g.fig.set_figwidth(5)
        g.fig.set_figheight(7)
        plt.title("Dispersion de la 1ere vp de "+findhist + "NbOccur="+ str(len(dpol["pol"+findhist])))
        g=sns.displot(alleig[findhist][0][1],label='Valeur propre 1',bins=len(str(len(alleig[name_type_final[q]][0][1])))*20,kde=True,color='y')
        g.fig.set_figwidth(5)
        g.fig.set_figheight(7)
        plt.title("Dispersion de la 2eme vp de "+findhist + "NbOccur="+ str(len(dpol["pol"+findhist])))
        g=sns.displot(alleig[findhist][0][2],label='Valeur propre 1',bins=len(str(len(alleig[name_type_final[q]][0][2])))*20,kde=True,color='g')
        g.fig.set_figwidth(5)
        g.fig.set_figheight(7)
        plt.title("Dispersion de la 3eme vp de "+findhist + "NbOccur="+ str(len(dpol["pol"+findhist])))
        plt.show()


    #création des histogramme de la dispersion des 3 vp pour chaque type
    #os.chdir(image_path)
    ##Tracé des histogrammes

    if(create_distr_eig==0):
        for q in range(len(name_type_final)):

            sns.set_style("whitegrid")
            g=sns.displot(data=alleig[name_type_final[q]][0][0],bins=len(str(len(alleig[name_type_final[q]][0][0])))*20,kde=True,color='red')
            g.fig.set_figwidth(5)
            g.fig.set_figheight(7)
            plt.title("Dispersion de la 1ere vp de "+name_type_final[q]+"NbOccur="+str(len(dpol["pol"+name_type_final[q]])))
            plt.savefig(name_type_final[q]+"_VP1_DISP_S66"+'.png')
            g=sns.displot(data=alleig[name_type_final[q]][0][1],bins=len(str(len(alleig[name_type_final[q]][0][1])))*20,kde=True,color='y')
            g.fig.set_figwidth(5)
            g.fig.set_figheight(7)
            plt.title("Dispersion de la 2eme vp de "+name_type_final[q]+"NbOccur="+str(len(dpol["pol"+name_type_final[q]])))
            plt.savefig(name_type_final[q]+"_VP2_DISP_S66"+'.png')
            g=sns.displot(data=alleig[name_type_final[q]][0][2],bins=len(str(len(alleig[name_type_final[q]][0][2])))*20,kde=True,color='g')
            g.fig.set_figwidth(5)
            g.fig.set_figheight(7)
            plt.title("Dispersion de la 3eme vp de "+name_type_final[q]+"NbOccur="+str(len(dpol["pol"+name_type_final[q]])))
            plt.savefig(name_type_final[q]+"_VP3_DISP_S66"+'.png')

        for i in range(len(d)):
            os.chdir(hist_path)
            os.makedirs(name_type_final[i]+"_DIST_EIG_S66")

        os.chdir(image_path)
        image=os.listdir()
        path=hist_path

        for x in range(len(image)):
            for i in range(len(name_type_final)):
                if image[x].startswith(name_type_final[i]):
                    shutil.copy(image[x], os.path.join(path,name_type_final[i]+"_DIST_EIG_S66"))


    ##différentes bibliothèque de parametres statistique
    #ecart type des eigval
    for i in range(len(name_type_final)):
        std[name_type_final[i]].append(np.std(alleig[name_type_final[i]][0][0]))
        std[name_type_final[i]].append(np.std(alleig[name_type_final[i]][0][1]))
        std[name_type_final[i]].append(np.std(alleig[name_type_final[i]][0][2]))


    tr = {}
    for i in range(len(name_type_final)):
        tr['trace'+name_type_final[i]] = []


    for i in range(len(eigv_pol)):
        for j in range(len(eigv_pol["eigv_pol"+name_type_final[i]])):
            b=((eigv_pol["eigv_pol"+name_type_final[i]][j][0][0]+eigv_pol["eigv_pol"+name_type_final[i]][j][0][1]+eigv_pol["eigv_pol"+name_type_final[i]][j][0][2])/3)
            tr['trace'+name_type_final[i]].append(b)



    #fonction qui permet de rentrer un intervalle numerique dune valeur porpre sur un type donné et renvoie les dossier ou on la trouve

    def eigv_to_atom(eigsup,eiginf,num_eigv,type_atm):
        for i in range(3,len(ligne2)):
            t=np.array([[ligne2[i][7],ligne2[i][10],ligne2[i][11]],[ligne2[i][10],ligne2[i][8],ligne2[i][12]],[ligne2[i][11],ligne2[i][12],ligne2[i][9]]],dtype=float)
            p=[ligne2[i][14],ligne2[i][15],ligne2[i][16],ligne2[i][17],ligne2[i][18],ligne2[i][19],ligne2[i][20],ligne2[i][21],lettre(ligne2[i][1]),ligne2[i][4],lettre(ligne2[i][5]),lettre(ligne2[i][6])]
            for j in range(len(d["type"+type_atm])):
                if (p == d["type"+type_atm][j] and (np.linalg.eigh(t)[0][num_eigv-1]>=eiginf and np.linalg.eigh(t)[0][num_eigv-1] <= eigsup)):
                    print('fichier',ligne2[i][31])
                    print('Numero ',ligne2[i][0])

    #Création d'une matrice de moyennes des éléments pour chaque type
    moy = {}
    for i in range(len(name_type_final)):
        moy[name_type_final[i]] = []

    for i in range(len(dpol)):
        for j in range(len(dpol["pol" + name_type_final[i]])):
            if j==0:
                axx=float(dpol["pol" + name_type_final[i]][j][0])
                ayy=float(dpol["pol" + name_type_final[i]][j][1])
                azz=float(dpol["pol" + name_type_final[i]][j][2])
                axy=float(dpol["pol" + name_type_final[i]][j][3])
                axz=float(dpol["pol" + name_type_final[i]][j][4])
                ayz=float(dpol["pol" + name_type_final[i]][j][5])
            else:
                axx=axx+float(dpol["pol" + name_type_final[i]][j][0])
                ayy=ayy+float(dpol["pol" + name_type_final[i]][j][1])
                azz=azz+float(dpol["pol" + name_type_final[i]][j][2])
                axy=axy+float(dpol["pol" + name_type_final[i]][j][3])
                axz=axz+float(dpol["pol" + name_type_final[i]][j][4])
                ayz=ayz+float(dpol["pol" + name_type_final[i]][j][5])

        axx=axx/(len(dpol["pol" + name_type_final[i]]))
        ayy=ayy/(len(dpol["pol" + name_type_final[i]]))
        azz=azz/(len(dpol["pol" + name_type_final[i]]))
        axy=axy/(len(dpol["pol" + name_type_final[i]]))
        axz=axz/(len(dpol["pol" + name_type_final[i]]))
        ayz=ayz/(len(dpol["pol" + name_type_final[i]]))
        ayx=axy
        azx=axz
        azy=ayz

        m=np.array([axx,axy,axz,ayx,ayy,azy,azx,azy,azz],dtype=float)
        mm=(m.reshape(3,3))
        moy[name_type_final[i]].append(mm)


    #Ecart type pour éléments des matrices pour chaque type

    ec = {}
    for i in range(len(name_type_final)):
        ec[name_type_final[i]] = []
    alpha_xx=[]
    alpha_yy=[]
    alpha_zz=[]
    alpha_xy=[]
    alpha_xz=[]
    alpha_yz=[]
    for i in range(len(dpol)):
        for j in range(len(dpol["pol" + name_type_final[i]])):
                axx=float(dpol["pol" + name_type_final[i]][j][0])
                alpha_xx.append(axx)
                ayy=float(dpol["pol" + name_type_final[i]][j][1])
                alpha_yy.append(ayy)
                azz=float(dpol["pol" + name_type_final[i]][j][2])
                alpha_zz.append(azz)
                axy=float(dpol["pol" + name_type_final[i]][j][3])
                alpha_xy.append(axy)
                axz=float(dpol["pol" + name_type_final[i]][j][4])
                alpha_xz.append(axz)
                ayz=float(dpol["pol" + name_type_final[i]][j][5])
                alpha_yz.append(ayz)
        alpha_xx_ec=np.std(alpha_xx)
        alpha_yy_ec=np.std(alpha_yy)
        alpha_zz_ec=np.std(alpha_zz)
        alpha_xy_ec=np.std(alpha_xy)
        alpha_xz_ec=np.std(alpha_xz)
        alpha_yz_ec=np.std(alpha_yz)
        alpha_yx_ec=np.std(alpha_xy)
        alpha_zx_ec=np.std(alpha_xz)
        alpha_zy_ec=np.std(alpha_yz)

        a=np.array([alpha_xx_ec,alpha_xy_ec,alpha_xz_ec,alpha_yx_ec,alpha_yy_ec,alpha_yz_ec,alpha_zx_ec,alpha_zy_ec,alpha_zz_ec])
        e_c=(a.reshape(3,3))
        ec[name_type_final[i]].append(e_c)

    def lettre(X):
        return str(re.sub('[0123456789_]','',X))


    ##test K2 AGOSTINO sur les vp( grandes occurences )
    #dictionnaire resultats du test
    agostinok2_eig={}
    for i in range(len(name_type_final)):
        agostinok2_eig[name_type_final[i]]=[]

    for i in range(len(name_type_final)):
        if (len(alleig[name_type_final[i]][0][0])>20):
            alpha=0.005
            a=scipy.stats.normaltest(alleig[name_type_final[i]][0][0])
            b=scipy.stats.normaltest(alleig[name_type_final[i]][0][1])
            c=scipy.stats.normaltest(alleig[name_type_final[i]][0][2])
            if a[1]>alpha:
                agostinok2_eig[name_type_final[i]].append('p='+str(a[1])+'  eig datta follows normal')
            else:
                agostinok2_eig[name_type_final[i]].append('Eig datta does not follow normal with alpha=5% ')
            if b[1]>alpha:
                agostinok2_eig[name_type_final[i]].append('p='+str(b[1])+'   eig datta follows normal')
            else:
                agostinok2_eig[name_type_final[i]].append('Eig datta does not follow normal with alpha=5%')
            if c[1]>alpha:
                agostinok2_eig[name_type_final[i]].append('p='+str(c[1])+' eig datta follows normal')
            else:
                agostinok2_eig[name_type_final[i]].append('Eig datta does not follow normal with alpha=5%')
        else:
            agostinok2_eig[name_type_final[i]].append('Nb eig for' + name_type_final[i]+' <20 : K2 test irelevant')


    ##création du POL_OUTPUT_TRANSFORM+TYPE avec le type atomique

    #ajout du type a ligne2
    for i in range(3,len(ligne2)):
        p=[ligne2[i][14],ligne2[i][15],ligne2[i][16],ligne2[i][17],ligne2[i][18],ligne2[i][19],ligne2[i][20],ligne2[i][21],lettre(ligne2[i][1]),ligne2[i][4],lettre(ligne2[i][5]),lettre(ligne2[i][6])]
        for k in range(len(d)):
            for j in range(len(d["type"+name_type_final[k]])):
                if (p==d["type"+name_type_final[k]][j]):
                    ligne2[i].append(name_type_final[k])

    os.chdir(path_to_txtfiles)
    with open("pol_OUTPUT_TRANSFORM+TYPE.txt",'w') as f:
        for i in range(3,len(ligne2)):
            f.write(ligne2[i][0]+'\t'+ligne2[i][0]+'\t'+ligne2[i][1]+'\t'+ligne2[i][2]+'\t'+ligne2[i][3]+'\t'+ligne2[i][4]+'\t'+ligne2[i][5]+'\t'+ligne2[i][6]+'\t'+ligne2[i][7]+'\t'+ligne2[i][8]+'\t'+ligne2[i][9]+'\t'+ligne2[i][10]+'\t'+ligne2[i][11]+'\t'+ligne2[i][12]+'\t'+ligne2[i][13]+'\t'+ligne2[i][14]+'\t'+ligne2[i][15]+'\t'+ligne2[i][16]+'\t'+ligne2[i][17]+'\t'+ligne2[i][18]+'\t'+ligne2[i][19]+'\t'+ligne2[i][20]+'\t'+ligne2[i][21]+'\t'+ligne2[i][22]+'\t'+ligne2[i][23]+'\t'+ligne2[i][24]+'\t'+ligne2[i][25]+'\t'+ligne2[i][26]+'\t'+ligne2[i][27]+'\t'+ligne2[i][28]+'\t'+ligne2[i][29]+'\t'+ligne2[i][30]+'\t'+ligne2[i][31]+'\t'+ligne2[i][32]+'\n')


    ##calcul des polarisabilité moléculaire
    os.chdir(path_to_the_bigset)
    #recupere les types atomiques présent dans un fichier
    typefile={}
    for i in range(len(fichiers)):
        typefile[fichiers[i]]=[]
    for i in range(3,len(ligne2)):
        for j in range(len(fichiers)):
            if (ligne2[i][31]==fichiers[j]):
                typefile[fichiers[j]].append(ligne2[i][32])
    #creation d'un dictionnaire qui regroupe toute les infos du pol output par fichiers
    file_infos={}
    for i in range(len(fichiers)):
        file_infos[fichiers[i]]=[]
    #remplissage du dictionnaire
    for i in range(3,len(ligne2)):
        for j in range(len(fichiers)):
            if ligne2[i][31]==fichiers[j]:
                file_infos[fichiers[j]].append(ligne2[i])
    #Diagonalisation des matrices moyennes elem par elem

    moy_vec_propre = {}
    test={}
    polar_test={}
    eig_pol_test={}
    tr_test={}
    for i in range(len(name_type_final)):
        moy_vec_propre[name_type_final[i]] = []
    for i in range(len(file_infos)):
        test[fichiers[i]]=[]
        polar_test[fichiers[i]]=[]
        eig_pol_test[fichiers[i]]=[]
        tr_test[fichiers[i]]=[]

    for i in range(len(moy)):
        a=np.array(moy[name_type_final[i]])
        b=np.linalg.eigh(a.reshape(3,3))
        moy_vec_propre[name_type_final[i]].append(b[1])
    ##Calcul des poalrisabilités moléculaires transférées

    moy_mol={}
    eig_moy_mol={}
    polar_mol={}
    tr_pol_mol={}
    delta_alpha={}
    delta_alpha_test={}
    for i in range(len(fichiers)):
        moy_mol[fichiers[i]]=[]
        eig_moy_mol[fichiers[i]]=[]
        polar_mol[fichiers[i]]=[]
        tr_pol_mol[fichiers[i]]=[]
        delta_alpha[fichiers[i]]=[]
        delta_alpha_test[fichiers[i]]=[]

    for i in range(len(file_infos)):
        for j in range(len(file_infos[fichiers[i]])):
            a=np.array([file_infos[fichiers[i]][j][22],file_infos[fichiers[i]][j][23],file_infos[fichiers[i]][j][24],file_infos[fichiers[i]][j][25],file_infos[fichiers[i]][j][26],file_infos[fichiers[i]][j][27],file_infos[fichiers[i]][j][28],file_infos[fichiers[i]][j][29],file_infos[fichiers[i]][j][30]], dtype=float)
            b=(a.reshape(3,3))
            base=np.transpose(b)
            tpol_mol=(np.dot(base,np.dot(moy[file_infos[fichiers[i]][j][32]][0],np.linalg.inv(base))))
            moy_mol[fichiers[i]].append(tpol_mol)
        s=np.zeros((3,3))
        for k in range(len(moy_mol[fichiers[i]])):
            s=np.add(s,moy_mol[fichiers[i]][k])
        polar_mol[fichiers[i]]=s
        eig_moy_mol[fichiers[i]].append(np.linalg.eigh(polar_mol[fichiers[i]])[0])
        tr_pol_mol[fichiers[i]]=((eig_moy_mol[fichiers[i]][0][0]+eig_moy_mol[fichiers[i]][0][1]+eig_moy_mol[fichiers[i]][0][2])/3)
        delta_alpha[fichiers[i]]=sqrt((3*np.trace(np.dot(polar_mol[fichiers[i]],polar_mol[fichiers[i]]))-(np.trace(polar_mol[fichiers[i]]))**2)/2)

    ##Calcul des poalrisabilités moléculaires théoriques

    polar_mol_th={}
    eig_th_mol={}
    tr_pol_mol_th={}
    th_mol={}
    delta_alpha_th={}

    for i in range(len(file_infos)):
        th_mol[fichiers[i]]=[]
        polar_mol_th[fichiers[i]]=[]
        eig_th_mol[fichiers[i]]=[]
        tr_pol_mol_th[fichiers[i]]=[]
        delta_alpha_th[fichiers[i]]=[]
        for j in range(len(file_infos[fichiers[i]])):
            a=np.array([file_infos[fichiers[i]][j][22],file_infos[fichiers[i]][j][23],file_infos[fichiers[i]][j][24],file_infos[fichiers[i]][j][25],file_infos[fichiers[i]][j][26],file_infos[fichiers[i]][j][27],file_infos[fichiers[i]][j][28],file_infos[fichiers[i]][j][29],file_infos[fichiers[i]][j][30]], dtype=float)
            b=(a.reshape(3,3))
            base=np.transpose(b)
            c=np.array([file_infos[fichiers[i]][j][7],file_infos[fichiers[i]][j][10],file_infos[fichiers[i]][j][11],file_infos[fichiers[i]][j][10],file_infos[fichiers[i]][j][8],file_infos[fichiers[i]][j][12],file_infos[fichiers[i]][j][11],file_infos[fichiers[i]][j][12],file_infos[fichiers[i]][j][9]], dtype=float)
            tenseur=(c.reshape(3,3))
            tenseur_th=(np.dot(base,np.dot(tenseur,np.linalg.inv(base))))
            th_mol[fichiers[i]].append(tenseur_th)
        s=np.zeros((3,3))
        for k in range(len(th_mol[fichiers[i]])):
            s=np.add(s,th_mol[fichiers[i]][k])
        polar_mol_th[fichiers[i]]=s
        eig_th_mol[fichiers[i]].append(np.linalg.eigh(polar_mol_th[fichiers[i]])[0])
        tr_pol_mol_th[fichiers[i]]=((eig_th_mol[fichiers[i]][0][0]+eig_th_mol[fichiers[i]][0][1]+eig_th_mol[fichiers[i]][0][2])/3)
        delta_alpha_th[fichiers[i]]=sqrt((3*np.trace(np.dot(polar_mol_th[fichiers[i]],polar_mol_th[fichiers[i]]))-(np.trace(polar_mol_th[fichiers[i]]))**2)/2)
#tarcé
    os.chdir(image_path)
    plt.figure()
    X=[]
    Y=[]
    a=0
    b=0
    for i in range(len(tr_pol_mol_th)):
        X.append(tr_pol_mol_th[fichiers[i]])
        Y.append(tr_pol_mol[fichiers[i]])
    a,b=np.polyfit(X,Y,1)
    coef = np.polyfit(X,Y,1)
    poly1d_fn = np.poly1d(coef)
    plt.title("Comparaison Trans/Theo pour le set de donnée")
    plt.xlabel('alpha theo')
    plt.ylabel('alpha trans')
    plt.text(0.15*(max(X)),0.75*max(poly1d_fn(X)),s=('y='+str(a.round(3))+'x+'+str(b.round(3))))
    plt.text(0.15*(max(X)),0.65*max(poly1d_fn(X)),s=('R²='+str(r2_score(X,Y).round(3))))
    plt.plot(X,Y,'r+', X, poly1d_fn(X),'-k')
    plt.savefig('alpha_.png')
    print('alpha compare saved in', image_path)
    print(' ')


    #évaluation de l'anisotropie des polarisabilités moléculaires
    os.chdir(image_path)
    X=[]
    Y=[]
    a=0
    b=0
    plt.figure()
    for i in range(len(delta_alpha)):
        Y.append(delta_alpha[fichiers[i]])
        X.append(delta_alpha_th[fichiers[i]])
    a,b=np.polyfit(X,Y,1)
    coef = np.polyfit(X,Y,1)
    poly1d_fn = np.poly1d(coef)
    plt.title("Anisotropie Trans/Theo pour le set de donnée")
    plt.xlabel('delta alpha theo')
    plt.ylabel('delta alpha trans')
    plt.xlim(0,max(X))
    plt.ylim(0,max(poly1d_fn(X)))
    plt.text(0.15*(max(X)),0.75*max(poly1d_fn(X)),s=('y='+str(a.round(3))+'x+'+str(b.round(3))))
    plt.text(0.15*(max(X)),0.65*max(poly1d_fn(X)),s=('R²='+str(r2_score(X,Y).round(3))))
    plt.plot(X,Y,'y+', X, poly1d_fn(X),'-k')
    plt.savefig('delta_alpha_.png')
    print('delta alpha compare saved in', image_path)
    print(' ')
    print(' ###  eigen value distribution tool callable with : find_distr(atomic_type_name) ###')
    print(' ')
    print('run complete in ',(time.time() - start_time))

    # ## Methode avec les vecteurs propres(!!!uncomplete)
    # #eigvect mean
    # eigvect_pol_mean={}
    # ortho_eigvect_pol_mean={}
    # for i in range(len(name_type_final)):
    #     eigvect_pol_mean[name_type_final[i]] = []
    #     ortho_eigvect_pol_mean[name_type_final[i]] = []
    #
    #
    # #fonction pour la moyenne dun vecteur propre (1,2,3) par type
    # def eigvmean(a,num_vect,b):
    #     vref=b[num_vect-1]
    #     som=[0]*3
    #     for i in range(0,len(a)):
    #         if np.vdot(vref,a[i][1][num_vect-1])>=0:
    #             som=som+a[i][1][num_vect-1]
    #         else:
    #             som=som-a[i][1][num_vect-1]
    #     return np.array(som)/len(a)
    #
    # for i in range(len(name_type_final)):
    #     eigvect_pol_mean[name_type_final[i]].append([eigvmean(eigv_pol["eigv_pol"+name_type_final[i]],1,moy_vec_propre[name_type_final[i]][0]),eigvmean(eigv_pol["eigv_pol"+name_type_final[i]],2,moy_vec_propre[name_type_final[i]][0]),eigvmean(eigv_pol["eigv_pol"+name_type_final[i]],3,moy_vec_propre[name_type_final[i]][0])])
    #
    #
    # #algo dorthogonalisation des vect propre (source : theo phd abridged, uncomplète )
    # def ortho_mean_eigvect(v1,v2):
    #     epsilon=1
    #     while epsilon>0.00002:
    #         v2_orth=np.array(v2-0.5*(np.vdot(v1,v2))*v1)
    #         v1_orth=np.array(v1-0.5*(np.vdot(v1,v2))*v2)
    #         e_1=v1_orth/(np.linalg.norm(v1_orth))
    #         e_2=v2_orth/(np.linalg.norm(v2_orth))
    #         epsilon = np.vdot(e_1,e_2)
    #         v1=v1_orth
    #         v2=v2_orth
    #     e_3=np.cross(e_1,e_2)/(np.linalg.norm(np.cross(e_1,e_2)))
    #     return [e_1,e_2,e_3]
    #
    #
    # for i in range(len(eigvect_pol_mean)):
    #     ortho_eigvect_pol_mean[name_type_final[i]].append(ortho_mean_eigvect(eigvect_pol_mean[name_type_final[i]][0][0],eigvect_pol_mean[name_type_final[i]][0][1]))
    #
    # #fonction moyenne des vp
    # def eigmean(a,num_vp):
    #     long=len(a)
    #     som=0
    #     for i in range(long):
    #         som=som+a[i][0][num_vp-1]
    #     return som/long
    #
    #
    # eigv_pol_mean={}
    # for i in range(len(name_type_final)):
    #     eigv_pol_mean[name_type_final[i]] = []
    #
    # for i in range(len(name_type_final)):
    #         eigv_pol_mean[name_type_final[i]].append([eigmean(eigv_pol["eigv_pol"+name_type_final[i]],1),eigmean(eigv_pol["eigv_pol"+name_type_final[i]],2),eigmean(eigv_pol["eigv_pol"+name_type_final[i]],3)])



    #
    # ##Sortie du fichier regroupant les infos sur les types atomiques
    # os.chdir(path_to_txtfiles)
    # with open("ATOM_TYPE_DATTA.txt",'w') as f:
    #     for i in range(len(name_type_final)):
    #             f.write('TYPE '+'\t'+name_type_final[i]+'\t'+str(len(dpol["pol"+name_type_final[i]]))+'\n')
    #             f.write('EIGEN VALUES MEAN'+'\n')
    #             f.write(str(eigv_pol_mean[name_type_final[i]][0][0])+'\t'+str(eigv_pol_mean[name_type_final[i]][0][1])+'\t'+str(eigv_pol_mean[name_type_final[i]][0][2])+'\n')
    #             f.write('AGOSTINO K2 TEST FOR THE EIG VALUES SERIES'+'\n')
    #             f.write(str(agostinok2_eig[name_type_final[i]][0])+'\n'+str(agostinok2_eig[name_type_final[i]][0])+'\n'+str(agostinok2_eig[name_type_final[i]][0])+'\n')
    #             f.write('STANDARD DEVIATION OF EIGEN VALUES'+'\n')
    #             f.write(str(std[name_type_final[i]][0])+'\t'+str(std[name_type_final[i]][1])+'\t'+str(std[name_type_final[i]][2])+'\n')
    #             f.write('ORTHO EIGEN VECTORS MEAN'+'\n')
    #             f.write(str(ortho_eigvect_pol_mean[name_type_final[i]][0][0][0])+'\t'+str(ortho_eigvect_pol_mean[name_type_final[i]][0][0][1])+'\t'+str(ortho_eigvect_pol_mean[name_type_final[i]][0][0][2])+'\n')
    #             f.write(str(ortho_eigvect_pol_mean[name_type_final[i]][0][1][0])+'\t'+str(ortho_eigvect_pol_mean[name_type_final[i]][0][1][1])+'\t'+str(ortho_eigvect_pol_mean[name_type_final[i]][0][1][2])+'\n')
    #             f.write(str(ortho_eigvect_pol_mean[name_type_final[i]][0][2][0])+'\t'+str(ortho_eigvect_pol_mean[name_type_final[i]][0][2][1])+'\t'+str(ortho_eigvect_pol_mean[name_type_final[i]][0][2][2])+'\n')
    #             f.write('TPOL MEAN ELEMENT BY ELEMENT'+'\n')
    #             f.write(str(moy[name_type_final[i]][0][0][0])+'\t'+str(moy[name_type_final[i]][0][0][1])+'\t'+str(moy[name_type_final[i]][0][0][2])+'\n')
    #             f.write(str(moy[name_type_final[i]][0][1][0])+'\t'+str(moy[name_type_final[i]][0][1][1])+'\t'+str(moy[name_type_final[i]][0][1][2])+'\n')
    #             f.write(str(moy[name_type_final[i]][0][2][0])+'\t'+str(moy[name_type_final[i]][0][2][1])+'\t'+str(moy[name_type_final[i]][0][2][2])+'\n')
    #             f.write('STANDARD DEVIATION TPOL MEAN ELEMENT'+'\n')
    #             f.write(str(ec[name_type_final[i]][0][0][0])+'\t'+str(ec[name_type_final[i]][0][0][1])+'\t'+str(ec[name_type_final[i]][0][0][2])+'\n')
    #             f.write(str(ec[name_type_final[i]][0][1][0])+'\t'+str(ec[name_type_final[i]][0][1][1])+'\t'+str(ec[name_type_final[i]][0][1][2])+'\n')
    #             f.write(str(ec[name_type_final[i]][0][2][0])+'\t'+str(ec[name_type_final[i]][0][2][1])+'\t'+str(ec[name_type_final[i]][0][2][2])+'\n')
    #             f.write('\n')
    ###########Sortie des fichiers ELMAM pour S66
    # os.chdir(path_to_txtfiles)
    # with open("Sortie_TPol.txt",'w') as f:
    #     for i in range(len(name_type_final)):
    #         f.write(name_type_final[i]+'\n')
    #         f.write(str(moy[name_type_final[i]][0][0][0].round(5)).format(num)+' '+str(moy[name_type_final[i]][0][1][1].round(5)).format(num)+' '+str(moy[name_type_final[i]][0][2][2].round(5)).format(num)+' '+str(moy[name_type_final[i]][0][0][1].round(5)).format(num)+' '+str(moy[name_type_final[i]][0][0][2].round(5)).format(num)+' '+str(moy[name_type_final[i]][0][1][2].round(5)).format(num)+'\n')
    #         f.write('\n')


    ligne8=[]
    os.chdir(elmam_path)
    with open ("pELMAM3_BASE.txt", "r") as file:
        data_ana=file.readlines()
        for li in data_ana :
            s = li.strip ("\n""\r""\t")
            l = s.split()
            ligne8.append (l)
    for i in range(len(ligne8)):
        for j in range(len(name_type_final)):
            if len(ligne8[i])!=0 and ligne8[i][0]=='ATOM'and ligne8[i][1]==name_type_final[j]:
                for k in range(20):
                    if 'ALPHA'in ligne8[i+k]:
                        z = '{:<07}'
                        data_ana[i+k]=('ALPHA'+'       '+z.format(moy[name_type_final[j]][0][0][0].round(5))+'  '+z.format(moy[name_type_final[j]][0][1][1].round(5))+'  '+z.format(moy[name_type_final[j]][0][2][2].round(5))+'  '+z.format(moy[name_type_final[j]][0][0][1].round(5))+'  '+z.format(moy[name_type_final[j]][0][0][2].round(5))+'  '+z.format(moy[name_type_final[j]][0][1][2].round(5))+'\n')

    with open ("S66_pELMAM3_BASE.txt", "w") as file2:
        file2.writelines(data_ana)




    #cas avec vecteurs propres (méthode uncomplète)
    # for i in range(len(file_infos)):
    #     for j in range(len(file_infos[fichiers[i]])):
    #         a=np.array([file_infos[fichiers[i]][j][22],file_infos[fichiers[i]][j][23],file_infos[fichiers[i]][j][24],file_infos[fichiers[i]][j][25],file_infos[fichiers[i]][j][26],file_infos[fichiers[i]][j][27],file_infos[fichiers[i]][j][28],file_infos[fichiers[i]][j][29],file_infos[fichiers[i]][j][30]], dtype=float)
    #         b=(a.reshape(3,3))
    #         c=np.array([file_infos[fichiers[i]][j][7],file_infos[fichiers[i]][j][10],file_infos[fichiers[i]][j][11],file_infos[fichiers[i]][j][10],file_infos[fichiers[i]][j][8],file_infos[fichiers[i]][j][12],file_infos[fichiers[i]][j][11],file_infos[fichiers[i]][j][12],file_infos[fichiers[i]][j][9]], dtype=float)
    #         tenseur=(c.reshape(3,3))
    #         diag_eig=np.array([eigv_pol_mean[file_infos[fichiers[i]][j][32]][0][0],0,0,0,eigv_pol_mean[file_infos[fichiers[i]][j][32]][0][1],0,0,0,eigv_pol_mean[file_infos[fichiers[i]][j][32]][0][2]],dtype=float)
    #         tpol=np.dot(np.linalg.inv(np.transpose(ortho_eigvect_pol_mean[file_infos[fichiers[i]][j][32]][0])),np.dot(np.reshape(diag_eig,(3,3)),np.transpose(ortho_eigvect_pol_mean[file_infos[fichiers[i]][j][32]][0])))
    #         T=np.dot(np.transpose(b),np.dot(tpol,np.linalg.inv(np.transpose(b))))
    #         test[fichiers[i]].append(T)
    #     s=np.zeros((3,3))
    #     for k in range(len(test[fichiers[i]])):
    #         s=np.add(s,test[fichiers[i]][k])
    #     polar_test[fichiers[i]]=s
    #     eig_pol_test[fichiers[i]].append(np.linalg.eigh(polar_test[fichiers[i]])[0])
    #     tr_test[fichiers[i]]=np.sum(eig_pol_test[fichiers[i]])/3
    #     delta_alpha_test[fichiers[i]]=sqrt((3*np.trace(np.dot(polar_test[fichiers[i]],polar_test[fichiers[i]]))-(np.trace(polar_test[fichiers[i]]))**2)/2)


    ##Transfer  du set d'Anna sur le data66
    Anna=['C301', 'C304', 'C305', 'C306', 'C307', 'C308T', 'C312', 'C314', 'C316', 'C4', 'C401', 'C402', 'C403', 'C407', 'C410', 'C411', 'C413T', 'C414', 'C408', 'CL-', 'H1 ', 'H101', 'H102', 'H105', 'H106', 'H107', 'H108', 'H103', 'H104', 'H109T', 'H110T', 'H111', 'H113', 'H111T', 'H112', 'H114T', 'H115', 'N201', 'N301', 'N302', 'N303', 'N304', 'N401', 'N402T', 'O102', 'O104', 'O105', 'O106', 'O1B', 'O201T', 'O202', 'O203', 'O204T', 'Oh1 ', 'P1', 'S202']

    data6=['C209', 'C301', 'C318', 'C306', 'C314', 'C315', 'C316', 'C401', 'C402T', 'C402', 'C405T', 'C406T', 'Cohhh', 'H101', 'H102', 'H106', 'H108', 'H103', 'H104', 'H111', 'H112', 'H114T', 'H115', 'H118', 'N201', 'N301', 'N302', 'N401', 'O102', 'O103T', 'O105', 'O201T', 'O202', 'O203']

    tyype=[]
    for i in range(len(Anna)):
        if Anna[i] not in tyype:
            tyype.append(Anna[i])

    for i in range(len(data6)):
        if data6[i] not in tyype:
            tyype.append(data6[i])


    common=set(Anna)&set(data6)
    #liste des types communs aux deux set de fichiers
    type_common=['N301', 'H114T', 'O201T', 'H106', 'N401', 'H108', 'C314', 'O105', 'C301', 'H101', 'N302', 'H115', 'C402', 'C316', 'C306', 'O102', 'H111', 'N201', 'O202', 'O203', 'H112', 'C401', 'H103', 'H104', 'H102']

    #liste des fichiers pouvant etre transférés avec cette liste de type atomiques communs aux deux sets
    allowed_file=[]
    not_allowed_file=[]
    a=[]

    for i in range(len(file_infos)):
        for j in range(len(file_infos[fichiers[i]])):
            if file_infos[fichiers[i]][j][32] not in type_common:
                a.append(fichiers[i])
    for i in range(len(a)):
        if a[i] not in not_allowed_file:
            not_allowed_file.append(a[i])
    for i in range(len(fichiers)):
        if fichiers[i] not in not_allowed_file:
            allowed_file.append(fichiers[i])

    # ##calcul des polarisabilités moléculaires sur les 20 allowed files
    # #!!! attention les varaibales ici ont été mise a zero, ne pas oublier de toggle comment cette partie apres l'avoir terminée
    # #calcul des 20 transférés
    # from ClassicoAnna import moy
    # moyy=moy
    # moy_mol={}
    # eig_moy_mol={}
    # polar_mol={}
    # tr_pol_mol={}
    # delta_alpha={}
    # delta_alpha_test={}
    # for i in range(len(allowed_file)):
    #     moy_mol[allowed_file[i]]=[]
    #     eig_moy_mol[allowed_file[i]]=[]
    #     polar_mol[allowed_file[i]]=[]
    #     tr_pol_mol[allowed_file[i]]=[]
    #     delta_alpha[allowed_file[i]]=[]
    #     delta_alpha_test[allowed_file[i]]=[]
    #
    # for i in range(len(allowed_file)):
    #     for j in range(len(file_infos[allowed_file[i]])):
    #         a=np.array([file_infos[allowed_file[i]][j][22],file_infos[allowed_file[i]][j][23],file_infos[allowed_file[i]][j][24],file_infos[allowed_file[i]][j][25],file_infos[allowed_file[i]][j][26],file_infos[allowed_file[i]][j][27],file_infos[allowed_file[i]][j][28],file_infos[allowed_file[i]][j][29],file_infos[allowed_file[i]][j][30]], dtype=float)
    #         b=(a.reshape(3,3))
    #         base=np.transpose(b)
    #         tpol_mol=(np.dot(base,np.dot(moyy[file_infos[allowed_file[i]][j][32]][0],np.linalg.inv(base))))
    #         moy_mol[allowed_file[i]].append(tpol_mol)
    #     s=np.zeros((3,3))
    #     for k in range(len(moy_mol[allowed_file[i]])):
    #         s=np.add(s,moy_mol[allowed_file[i]][k])
    #     polar_mol[allowed_file[i]]=s
    #     eig_moy_mol[allowed_file[i]].append(np.linalg.eigh(polar_mol[allowed_file[i]])[0])
    #     tr_pol_mol[allowed_file[i]]=((eig_moy_mol[allowed_file[i]][0][0]+eig_moy_mol[allowed_file[i]][0][1]+eig_moy_mol[allowed_file[i]][0][2])/3)
    #     delta_alpha[allowed_file[i]]=sqrt((3*np.trace(np.dot(polar_mol[allowed_file[i]],polar_mol[allowed_file[i]]))-(np.trace(polar_mol[allowed_file[i]]))**2)/2)
    #
    # #calcul des 20 thériques
    #
    # polar_mol_th={}
    # eig_th_mol={}
    # tr_pol_mol_th={}
    # th_mol={}
    # delta_alpha_th={}
    #
    # for i in range(len(allowed_file)):
    #     th_mol[allowed_file[i]]=[]
    #     polar_mol_th[allowed_file[i]]=[]
    #     eig_th_mol[allowed_file[i]]=[]
    #     tr_pol_mol_th[allowed_file[i]]=[]
    #     delta_alpha_th[allowed_file[i]]=[]
    #     for j in range(len(file_infos[allowed_file[i]])):
    #         a=np.array([file_infos[allowed_file[i]][j][22],file_infos[allowed_file[i]][j][23],file_infos[allowed_file[i]][j][24],file_infos[allowed_file[i]][j][25],file_infos[allowed_file[i]][j][26],file_infos[allowed_file[i]][j][27],file_infos[allowed_file[i]][j][28],file_infos[allowed_file[i]][j][29],file_infos[allowed_file[i]][j][30]], dtype=float)
    #         b=(a.reshape(3,3))
    #         base=np.transpose(b)
    #         c=np.array([file_infos[allowed_file[i]][j][7],file_infos[allowed_file[i]][j][10],file_infos[allowed_file[i]][j][11],file_infos[allowed_file[i]][j][10],file_infos[allowed_file[i]][j][8],file_infos[allowed_file[i]][j][12],file_infos[allowed_file[i]][j][11],file_infos[allowed_file[i]][j][12],file_infos[allowed_file[i]][j][9]], dtype=float)
    #         tenseur=(c.reshape(3,3))
    #         tenseur_th=(np.dot(base,np.dot(tenseur,np.linalg.inv(base))))
    #         th_mol[allowed_file[i]].append(tenseur_th)
    #     s=np.zeros((3,3))
    #     for k in range(len(th_mol[allowed_file[i]])):
    #         s=np.add(s,th_mol[allowed_file[i]][k])
    #     polar_mol_th[allowed_file[i]]=s
    #     eig_th_mol[allowed_file[i]].append(np.linalg.eigh(polar_mol_th[allowed_file[i]])[0])
    #     tr_pol_mol_th[allowed_file[i]]=((eig_th_mol[allowed_file[i]][0][0]+eig_th_mol[allowed_file[i]][0][1]+eig_th_mol[allowed_file[i]][0][2])/3)
    #     delta_alpha_th[allowed_file[i]]=sqrt((3*np.trace(np.dot(polar_mol_th[allowed_file[i]],polar_mol_th[allowed_file[i]]))-(np.trace(polar_mol_th[allowed_file[i]]))**2)/2)
    #
    # ##tracé rapide de la comparaison TRANSFERE/THEO + ANISOTROPIE: ELEM ELEM sur les 15 dimeres autorisé --------TRANS BY ANNA--------
    # os.chdir(r'C:\Users\alban\Desktop\Stages L3_M1\Stage M1\Python\Graph_Disp_Eig\ELEM_ELEM')
    # X=[]
    # Y=[]
    # a=0
    # b=0
    # for i in range(len(tr_pol_mol_th)):
    #     X.append(tr_pol_mol_th[allowed_file[i]])
    #     Y.append(tr_pol_mol[allowed_file[i]])
    # a,b=np.polyfit(X,Y,1)
    # coef = np.polyfit(X,Y,1)
    # poly1d_fn = np.poly1d(coef)
    # plt.title("Comparaison Trans/Theo pour le set66 transféré par ELMAMA_Anna")
    # plt.xlabel('alpha theo')
    # plt.ylabel('alpha trans')
    # plt.text(4,17,s=('y='+str(a.round(3))+'x+'+str(b.round(3))))
    # plt.text(4,16,s=('R²='+str(r2_score(X,Y).round(3))))
    # plt.plot(X,Y,'r+', X, poly1d_fn(X),'-k')
    # plt.savefig('alpha_S66Trans_ByAnna.png')
    # plt.show()
    # plt.close()
    # #évaluation de l'anisotropie des polarisabilités moléculaires
    # X=[]
    # Y=[]
    # a=0
    # b=0
    # for i in range(len(delta_alpha)):
    #     Y.append(delta_alpha[allowed_file[i]])
    #     X.append(delta_alpha_th[allowed_file[i]])
    # a,b=np.polyfit(X,Y,1)
    # coef = np.polyfit(X,Y,1)
    # poly1d_fn = np.poly1d(coef)
    # plt.title("Anisotropie Trans/Theo pour le set66 transféré par ELMAMA_Anna")
    # plt.xlabel('delta alpha theo')
    # plt.ylabel('delta alpha trans')
    # plt.text(2,11,s=('y='+str(a.round(3))+'x+'+str(b.round(3))))
    # plt.text(2,10,s=('R²='+str(r2_score(X,Y).round(3))))
    # plt.plot(X,Y,'y+', X, poly1d_fn(X),'-k')
    # plt.savefig('delta_alpha_S66Trans_ByAnna.png')
    # plt.show()
    # plt.close()

    ##tracé rapide de la comparaison TRANSFERE/THEO + ANISOTROPIE: ELEM ELEM sur les 15 dimeres autorisé -------- TRANS BY S66--------
    if toogle_dim==1:
        moy_mol={}
        eig_moy_mol={}
        polar_mol={}
        tr_pol_mol={}
        delta_alpha={}
        delta_alpha_test={}
        for i in range(len(allowed_file)):
            moy_mol[allowed_file[i]]=[]
            eig_moy_mol[allowed_file[i]]=[]
            polar_mol[allowed_file[i]]=[]
            tr_pol_mol[allowed_file[i]]=[]
            delta_alpha[allowed_file[i]]=[]
            delta_alpha_test[allowed_file[i]]=[]

        for i in range(len(allowed_file)):
            for j in range(len(file_infos[allowed_file[i]])):
                a=np.array([file_infos[allowed_file[i]][j][22],file_infos[allowed_file[i]][j][23],file_infos[allowed_file[i]][j][24],file_infos[allowed_file[i]][j][25],file_infos[allowed_file[i]][j][26],file_infos[allowed_file[i]][j][27],file_infos[allowed_file[i]][j][28],file_infos[allowed_file[i]][j][29],file_infos[allowed_file[i]][j][30]], dtype=float)
                b=(a.reshape(3,3))
                base=np.transpose(b)
                tpol_mol=(np.dot(base,np.dot(moy[file_infos[allowed_file[i]][j][32]][0],np.linalg.inv(base))))
                moy_mol[allowed_file[i]].append(tpol_mol)
            s=np.zeros((3,3))
            for k in range(len(moy_mol[allowed_file[i]])):
                s=np.add(s,moy_mol[allowed_file[i]][k])
            polar_mol[allowed_file[i]]=s
            eig_moy_mol[allowed_file[i]].append(np.linalg.eigh(polar_mol[allowed_file[i]])[0])
            tr_pol_mol[allowed_file[i]]=((eig_moy_mol[allowed_file[i]][0][0]+eig_moy_mol[allowed_file[i]][0][1]+eig_moy_mol[allowed_file[i]][0][2])/3)
            delta_alpha[allowed_file[i]]=sqrt((3*np.trace(np.dot(polar_mol[allowed_file[i]],polar_mol[allowed_file[i]]))-(np.trace(polar_mol[allowed_file[i]]))**2)/2)

        #Calcul des poalrisabilités moléculaires théoriques

        polar_mol_th={}
        eig_th_mol={}
        tr_pol_mol_th={}
        th_mol={}
        delta_alpha_th={}

        for i in range(len(allowed_file)):
            th_mol[allowed_file[i]]=[]
            polar_mol_th[allowed_file[i]]=[]
            eig_th_mol[allowed_file[i]]=[]
            tr_pol_mol_th[allowed_file[i]]=[]
            delta_alpha_th[allowed_file[i]]=[]
            for j in range(len(file_infos[allowed_file[i]])):
                a=np.array([file_infos[allowed_file[i]][j][22],file_infos[allowed_file[i]][j][23],file_infos[allowed_file[i]][j][24],file_infos[allowed_file[i]][j][25],file_infos[allowed_file[i]][j][26],file_infos[allowed_file[i]][j][27],file_infos[allowed_file[i]][j][28],file_infos[allowed_file[i]][j][29],file_infos[allowed_file[i]][j][30]], dtype=float)
                b=(a.reshape(3,3))
                base=np.transpose(b)
                c=np.array([file_infos[allowed_file[i]][j][7],file_infos[allowed_file[i]][j][10],file_infos[allowed_file[i]][j][11],file_infos[allowed_file[i]][j][10],file_infos[allowed_file[i]][j][8],file_infos[allowed_file[i]][j][12],file_infos[allowed_file[i]][j][11],file_infos[allowed_file[i]][j][12],file_infos[allowed_file[i]][j][9]], dtype=float)
                tenseur=(c.reshape(3,3))
                tenseur_th=(np.dot(base,np.dot(tenseur,np.linalg.inv(base))))
                th_mol[allowed_file[i]].append(tenseur_th)
            s=np.zeros((3,3))
            for k in range(len(th_mol[allowed_file[i]])):
                s=np.add(s,th_mol[allowed_file[i]][k])
            polar_mol_th[allowed_file[i]]=s
            eig_th_mol[allowed_file[i]].append(np.linalg.eigh(polar_mol_th[allowed_file[i]])[0])
            tr_pol_mol_th[allowed_file[i]]=((eig_th_mol[allowed_file[i]][0][0]+eig_th_mol[allowed_file[i]][0][1]+eig_th_mol[allowed_file[i]][0][2])/3)
            delta_alpha_th[allowed_file[i]]=sqrt((3*np.trace(np.dot(polar_mol_th[allowed_file[i]],polar_mol_th[allowed_file[i]]))-(np.trace(polar_mol_th[allowed_file[i]]))**2)/2)

        #tracé
        os.chdir(image_path)
        X=[]
        Y=[]
        a=0
        b=0
        for i in range(len(tr_pol_mol_th)):
            X.append(tr_pol_mol_th[allowed_file[i]])
            Y.append(tr_pol_mol[allowed_file[i]])
        a,b=np.polyfit(X,Y,1)
        coef = np.polyfit(X,Y,1)
        poly1d_fn = np.poly1d(coef)
        plt.title("Comparaison Trans/Theo pour le set66 transféré par ELMAMA_S66")
        plt.xlabel('alpha theo')
        plt.ylabel('alpha trans')
        plt.text(4,17,s=('y='+str(a.round(3))+'x+'+str(b.round(3))))
        plt.text(4,16,s=('R²='+str(r2_score(X,Y).round(3))))
        plt.plot(X,Y,'r+', X, poly1d_fn(X),'-k')
        plt.savefig('alpha_S66Trans_Bys66.png')
        plt.show()
        plt.close()
        #évaluation de l'anisotropie des polarisabilités moléculaires
        os.chdir(image_path)
        X=[]
        Y=[]
        a=0
        b=0
        for i in range(len(delta_alpha)):
            Y.append(delta_alpha[allowed_file[i]])
            X.append(delta_alpha_th[allowed_file[i]])
        a,b=np.polyfit(X,Y,1)
        coef = np.polyfit(X,Y,1)
        poly1d_fn = np.poly1d(coef)
        plt.title("Anisotropie Trans/Theo pour le set66 transféré par ELMAMA_S66")
        plt.xlabel('delta alpha theo')
        plt.ylabel('delta alpha trans')
        plt.text(2,11,s=('y='+str(a.round(3))+'x+'+str(b.round(3))))
        plt.text(2,10,s=('R²='+str(r2_score(X,Y).round(3))))
        plt.plot(X,Y,'y+', X, poly1d_fn(X),'-k')
        plt.savefig('delta_alpha_S66Trans_Bys66.png')
        plt.show()
        plt.close()












