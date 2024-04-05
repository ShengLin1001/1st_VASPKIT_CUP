# coding: utf-8

# Name: Mengnan Cui
# E-mail: mncui88@126.com
# Usage: This program is suitable for heterojunctions with the same lattice angle of two 2D materials and different lattice constants
# Version: 1.0
# Data: 2019.07.31



from pymatgen.core.surface import SlabGenerator, generate_all_slabs, Structure, Lattice, get_d
from pymatgen.io.vasp.inputs import Incar, Kpoints, Kpoints_supported_modes, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Dynmat, Outcar, Oszicar
# Import the neccesary tools for making a Wulff shape

import os
import shutil
import re
import time
import subprocess
import numpy as np
import pandas as pd
import math
from fractions import Fraction

np.set_printoptions(suppress=True)



def cons(poscar):
    x = np.array(poscar[2][:],float)
    y = np.array(poscar[2][:],float)
    z = np.array(poscar[2][:],float)
    a = math.sqrt(x[0]**2+x[1]**2+x[2]**2)
    b = math.sqrt(y[0]**2+y[1]**2+y[2]**2)
    c = math.sqrt(z[0]**2+z[1]**2+z[2]**2)
    return a
# poscar must splited
def natom(poscar):
    natom = sum([item for item in [float(i) for i in poscar[6][:]]])
    return natom
# get x,y of supercell a 
def mismatch(le,a1,a2,natom1,natom2):
    global mtype
    mtype = input("Which slab choosed as base(please input 'A' , 'B' or 'AVERAGE':  ")
    mat = []
    if mtype == 'A':
        for i in range(le):
            x = i + 1
            for h in range(le):
                y = h + 1
                mat.append(x)
                mat.append(y)
                mat.append(x*natom1+y*natom2)
                mat.append(x*a1)
                mat.append(x*a1)
                mat.append(abs(x*a1-y*a2)/(x*a1))  
    elif mtype == 'B':
        for i in range(le):
            x = i + 1
            for h in range(le):
                y = h + 1
                mat.append(x)
                mat.append(y)
                mat.append(x*natom1+y*natom2)
                mat.append(y*a2)
                mat.append(y*a2)
                mat.append(abs(x*a1-y*a2)/(y*a2))              
    else: 
        for i in range(le):
            x = i + 1
            for h in range(le):
                y = h + 1
                mat.append(x)
                mat.append(y)
                mat.append(x*natom1+y*natom2)
                mat.append((x*a1+y*a2)/2)
                mat.append((x*a1+y*a2)/2)
                mat.append(2*abs(x*a1-y*a2)/(x*a1+y*a2))
    return mat
# thickness of slab and max_coordination ,min_coordination of slabs
def thickness(poscar):
    
    direct = np.array([poscar[2],poscar[3],poscar[4]],float)
    natoms = int(natom(poscar))
    c = float(poscar[4][2])
    coor = []
    for i in range(int(natoms)):
        line = i + 8
        lines = np.array(poscar[line][0:3],float)
        cart = np.dot(direct.T,lines)
        coor.append(cart)
    coor = np.array(coor,float)
    coor = coor.reshape(natoms,3)
    coord = coor[coor[:,2].argsort()]
    max_coor = coord[natoms-1,:]
    min_coor = coord[0,:]
    thick = (max_coor[2]-min_coor[2])
    return min_coor, max_coor, thick, coor


# Matching slabs and output poscar
def match(poscar1, poscar2, mtype, vector):
    dist = input("Please input the distance between two layers: ")
    vacc = input("Please input the vaccuum thickness: ")
    pos_a = [line.strip(' ').split() for line in open(poscar1)]
    pos_b = [line.strip(' ').split() for line in open(poscar2)]

    direct_a = np.array([pos_a[2],pos_a[3],pos_a[4]],float)
    direct_b = np.array([pos_a[2],pos_a[3],pos_a[4]],float)

    min1, max1, thick1, coor1 = thickness(pos_a)
    min2, max2, thick2, coor2 = thickness(pos_b)
    lengthc = thick1 + thick2 + float(vacc) + float(dist)
    # scale constant of b
    if mtype == 'A':
        print(vector)
        vector = np.dot(direct_a.T, vector)
        print(vector)
        ratio = float(pos_a[2][0])/float(pos_b[2][0])
        direct_C = np.array([[0,0,lengthc]])
        file = open('Heterojunction_A.vasp','w')
        atomtype = pos_a[5] +pos_b[5]
        atomnum = pos_a[6]+pos_b[6]
        # write poscar
        file.write("Heterojunction \n 1.0 \n")
        for i in direct_a[0:2]:
            file.write("{:10f}\t{:10f}\t{:10f}\n".format(i[0],i[1],i[2]))
        for i in direct_C[0,:]:
            file.write("{:10f}\t".format(i))
        file.write("\n")
        for i in atomtype[:]:
            file.write(str(i)+' ')
        file.write("\n")
        for i in atomnum[:]:
            file.write(str(i)+' ')
        file.write("\n")
        file.write("Cartesian \n")
        for i in coor1[:]:
            for h in i[:]:
                file.write("{:10f}\t".format(h))
            file.write("\n")
        for i in coor2[:]:
            i = i + vector
            file.write("{:10f}\t{:10f}\t{:10f}\t".format(i[0]*ratio,i[1]*ratio,i[2]+max1[2]-min2[2]+float(dist)))
            file.write("\n")
        file.close()
        print('\n' + "POSCAR file of Heterojunction has been generated, named as Heterojunction_A.vasp")
    
    elif mtype == 'B':
        vector = np.dot(direct_b.T, vector)
        ratio = float(pos_b[2][0])/float(pos_a[2][0])
        direct_C = np.array([[0,0,lengthc]])

        file = open('Heterojunction_B.vasp','w')
        atomtype = pos_a[5] +pos_b[5]
        atomnum = pos_a[6]+pos_b[6]
        # write poscar
        file.write("Heterojunction \n 1.0 \n")
        for i in direct_b[0:2]:
            file.write("{:10f}\t{:10f}\t{:10f}\n".format(i[0],i[1],i[2]))
        for i in direct_C[0,:]:
            file.write("{:10f}\t".format(i))
        file.write("\n")
        for i in atomtype[:]:
            file.write(str(i)+' ')
        file.write("\n")
        for i in atomnum[:]:
            file.write(str(i)+' ')
        file.write("\n")
        file.write("Cartesian \n")
        for i in coor1[:]:
            file.write("{:10f}\t".format(i[0]*ratio))
            file.write("{:10f}\t".format(i[1]*ratio))
            file.write("{:10f}\t".format(i[2]))
            file.write("\n")        
        for i in coor2[:]:
            i = i + vector
            file.write("{:10f}\t{:10f}\t".format(i[0],i[1]))
            file.write("{:10f}".format(float(i[2])+max1[2]-min2[2]+float(dist)))
            file.write("\n")
        file.close()
        print('\n' + "POSCAR file of Heterojunction has been generated, named as Heterojunction_B.vasp")
    else:
        
        ratio_a = (float(pos_b[2][0])+float(pos_a[2][0]))/(2*float(pos_a[2][0]))
        ratio_b = (float(pos_b[2][0])+float(pos_a[2][0]))/(2*float(pos_b[2][0]))

        direct_AB = (direct_a[0:3,:] + direct_b[0:3,:])/2
        direct_C = np.array([0,0,lengthc])
        vector = np.dot(direct_AB.T, vector)

        file = open('Heterojunction_AVERAGE.vasp','w')
        atomtype = pos_a[5] +pos_b[5]
        atomnum = pos_a[6]+pos_b[6]
        # write poscar
        file.write("Heterogeneous \n 1.0 \n")
        for i in direct_AB[0:2]:
            file.write("{:10f}\t{:10f}\t{:10f}\n".format(i[0],i[1],i[2]))
        
        for i in direct_C[:]:
            file.write("{:10f}\t".format(i))
        file.write("\n")
    
        for i in atomtype[:]:
            file.write(str(i)+' ')
        file.write("\n")
    
        for i in atomnum[:]:
            file.write(str(i)+' ')
        file.write("\n")
        file.write("Cartesian \n")
    
        for i in coor1[:]:
            file.write("{:10f}\t".format(i[0]*ratio_a))
            file.write("{:10f}\t".format(i[1]*ratio_a))
            file.write("{:10f}\t".format(i[2]))
            file.write("\n")   
        
        for i in coor2[:]:
            i = i + vector
            file.write("{:10f}\t{:10f}\t".format(i[0]*ratio_b,i[1]*ratio_b))
            file.write("{:10f}".format(float(i[2])+max1[2]-min2[2]+float(dist)))
            file.write("\n")
        file.close()
        print('\n' + "POSCAR file of Heterojunction has been generated, named as Heterojunction_AVERAGE.vasp")


# Loading structure

structure1 = input("Please enter the file name of underlying material(eg.02_MoS2_POSCAR): ")
structure2 = input("Please enter the file name of upper material(eg.02_Graphene_POSCAR): ")
pos1 =[line.strip(' ').split() for line in open(structure1)]
pos2 =[line.strip(' ').split() for line in open(structure2)]
a1 = cons(pos1)
a2 = cons(pos2)
natom1 = natom(pos1)
natom2 = natom(pos2)
print('\n'+"Structure Information: ",'\n'
+str(structure1)+'\n'
+"Lattice constant: " + str(a1) + '\n'
+"Num of atom: " + str(natom1) + '\n'
+str(structure2)+'\n'
+"Lattice constant: " + str(a2) + '\n'
+"Num of atom: " + str(natom2) + '\n'
)



# mis-match
num_mis = float(input("Please enter the mis-match range(e.g. 0.05）:"))
n_list = 15
mat = np.array(mismatch(n_list,a1,a2,natom1,natom2))
#mat = np.around(mat, decimals = 2)
mat = mat.reshape(n_list**2,6)


mismat = mat[mat[:,5]<num_mis]
mismat = mismat[mismat[:,2].argsort()]
print("num \tNx \tNy \tN_atom \tLattice_axb \tmis-match")
h = 0
for i in mismat[:]:
    h = h + 1
    print("{}\t {} \t{} \t{} \t {:.2f} x {:.2f}\t{:.5f}".format(h,int(i[0]),int(i[1]),int(i[2]),i[3],i[4],i[5]))


# match cell
#lines = input("Which lines of mis-match data: ")

num = int(input("Please enter the number of line(e.g. 1): ")) -1
lines = mismat[num,:]
stru1 = Poscar.from_file(structure1).structure
stru2 = Poscar.from_file(structure2).structure
stru1.make_supercell([[lines[0],0,0],[0,lines[0],0],[0,0,1]])
stru2.make_supercell([[lines[1],0,0],[0,lines[1],0],[0,0,1]])
#slabs1 = SlabGenerator(stru1,[0,0,1],min_slab_size=5,min_vacuum_size=10,in_unit_planes=False).get_slabs({("Mo","S"):2.5})
#slabs2 = SlabGenerator(stru2,[0,0,1],min_slab_size=5,min_vacuum_size=10,in_unit_planes=False).get_slabs({("Mo","S"):2.5})
filename1 = "POSCAR1"
filename2 = "POSCAR2"

stru1.to('poscar',filename1)
stru2.to('poscar',filename2)


vector = np.array([Fraction(s) for s in input("Please enter the translation vector of the upper material, separated by space(e.g. -1/3 1/3 0): ").split()],float)
match(filename1, filename2,mtype,vector)

print('\n' +
"Name: Mengnan Cui" + '\n' +
"E-mail: mncui88@126.com" + '\n'+
"Version: 1.0" + '\n' +
"Usage: This program is suitable for heterojunctions with the same lattice angle of two 2D materials and different lattice constants" + '\n' +
"Data: 2019.07.31"  + '\n'
"Thanks for your attention")
