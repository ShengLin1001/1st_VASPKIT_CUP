#!/usr/bin/env python3
from numpy import *
import numpy as np
import math
import os
import copy

os.system('cp POSCAR11 POSCAR')
os.system('(echo 406; echo 1; echo 4061) | vaspkit | grep 123')
os.system('mv POSCAR_REV POSCAR11') 
os.system('cp POSCAR22 POSCAR')
os.system('(echo 406; echo 1; echo 4061) | vaspkit | grep 123')
os.system('mv POSCAR_REV POSCAR22') 
def read_first(POSCAR11_or_22): 
    fp1 = open(POSCAR11_or_22, mode='r')
    fp1.readline()
    fractor = float(fp1.readline())
    la = [[] for m in range(3)]
    for i in range(3):
        la[i][:] = list(fp1.readline().split())
        la[i][0:3] = list(map(float, la[i][0:3]))
    la[0][0] = fractor * la[0][0]
    la[0][1] = fractor * la[0][1]
    la[0][2] = fractor * la[0][2]
    la[1][0] = fractor * la[1][0]
    la[1][1] = fractor * la[1][1]
    la[1][2] = fractor * la[1][2]
    la[2][0] = fractor * la[2][0]
    la[2][1] = fractor * la[2][1]
    la[2][2] = fractor * la[2][2]
    C = la[2][2]
    ele = fp1.readline()
    ele = ele.strip('\n').split()
    atom_num = list(fp1.readline().split())
    atom_num = list(map(int, atom_num))
    fp1.readline()
    atos = fp1.readlines()
    row = len(atos)
    fp1.close()
    
    fp2 = open(POSCAR11_or_22, mode='r')
    atos_ele_C_sort = [[] for m in range(row)]
    atos_ele_C = [[] for m in range(row)]
    for i in range(8):
        fp2.readline()
    for i in range(row):
        atos_ele_C_sort[i][:] = list(fp2.readline().split())
        atos_ele_C_sort[i][0:3] = list(map(float, atos_ele_C_sort[i][0:3]))
    for i in range(row):
        atos_ele_C[i][:] = atos_ele_C_sort[i][:]
    atos_ele_C_sort.sort(key=lambda x: x[2], reverse=True)
    return fractor, la, ele, atom_num, atos_ele_C, C, atos_ele_C_sort 
    
    fp2.close()

[fractor_A, la_A, ele_A, atom_num_A, atos_ele_C_A, C_A, atos_ele_C_sort_A] = read_first("POSCAR11")
[fractor_B, la_B, ele_B, atom_num_B, atos_ele_C_B, C_B, atos_ele_C_sort_B] = read_first("POSCAR22")

def dividing_layer(atos_ele_C_sort_A_or_B): 
    z_coord = [0 for m in range(len(atos_ele_C_sort_A_or_B))]
    for i in range(len(atos_ele_C_sort_A_or_B)):
        z_coord[i] = atos_ele_C_sort_A_or_B[i][2]
    layer_coors = sorted(set(z_coord),key=z_coord.index) 
    k = 0
    layer = [[] for o1 in range(len(layer_coors))]
    for layer_coor in layer_coors:
        for i in range(len(atos_ele_C_sort_A_or_B)):
            if abs(atos_ele_C_sort_A_or_B[i][2] - layer_coor) < 0.01:   
                layer[k].append(atos_ele_C_sort_A_or_B[i][:])
        k += 1
    return layer
layer_A = dividing_layer(atos_ele_C_sort_A)
layer_B = dividing_layer(atos_ele_C_sort_B)
layer_A1 = copy.deepcopy(layer_A)
layer_B1 = copy.deepcopy(layer_B)

def cal_dis(coor1, coor2, q, w): 
    distance = ((coor1[q][0] - coor2[w][0]) ** 2 + (coor1[q][1] - coor2[w][1]) ** 2 + (coor1[q][2] - coor2[w][2]) ** 2) ** 0.5
    return distance
def cal_bond_layer(layer_A1_or_B):   
    bond_le_layer = [[] for i in range(len(layer_A1_or_B) - 1)]
    for i in range(len(layer_A1_or_B) - 1):
        for m in range(len(layer_A1_or_B[i])):
            for n in range(len(layer_A1_or_B[i + 1])):
                bond_le_layer[i].append(cal_dis(layer_A1_or_B[i], layer_A1_or_B[i + 1], m, n))
    bond_layer = []
    for i in range(len(bond_le_layer)):
        bond_layer.append(min(bond_le_layer[i]))
    return bond_layer


def cartesian_to_direct(layer_A1_or_B, lattice1, lattice2): 
    transit = [[] for i in range(len(layer_A1_or_B))]
    for i in range(len(layer_A1_or_B)):
        for m in range(len(layer_A1_or_B[i])):
            transit[i].append(layer_A1_or_B[i][m]) 
    lattice = [[] for i in range(len(lattice1))]
    for i in range(len(lattice1) - 1):
        lattice[i].append(lattice1[i][0])
        lattice[i].append(lattice1[i][1])
        lattice[i].append(lattice1[i][2])
    lattice[2].append(lattice1[2][0])
    lattice[2].append(lattice1[2][1])
    lattice[2].append(lattice1[2][2])
    lattice_tran_mat = np.mat(np.transpose(lattice))
    convmat = np.linalg.inv(lattice_tran_mat)
    latt_inv = convmat.tolist()
    layer_A_or_B_D = [[] for i in range(len(transit))]
    for i in range(len(transit)):
        for m in range(len(transit[i])):
            layer_A_or_B_D[i].append([(latt_inv[0][0] * transit[i][m][0] + latt_inv[0][1] * transit[i][m][1] + latt_inv[0][2] * transit[i][m][2]), \
            (latt_inv[1][0] * transit[i][m][0] + latt_inv[1][1] * transit[i][m][1] + latt_inv[1][2] * transit[i][m][2]), \
            (latt_inv[2][0] * transit[i][m][0] + latt_inv[2][1] * transit[i][m][1] + latt_inv[2][2] * transit[i][m][2]), \
            (layer_A1_or_B[i][m][3])])
    return layer_A_or_B_D

def direct_to_cartesian(layer_A_or_B_D,  lattice1,lattice2): 
    layer_AorB_C = copy.deepcopy(layer_A_or_B_D)
    for i in range(len(layer_A_or_B_D)):
        for m in range(len(layer_A_or_B_D[i])):
            layer_AorB_C[i][m][0] = (layer_A_or_B_D[i][m][0] * lattice1[0][0] + layer_A_or_B_D[i][m][1] \
                       * lattice1[1][0] + layer_A_or_B_D[i][m][2] * lattice2[2][0])
            layer_AorB_C[i][m][1] = (layer_A_or_B_D[i][m][0] * lattice1[0][1] + layer_A_or_B_D[i][m][1] \
                       * lattice1[1][1] + layer_A_or_B_D[i][m][2] * lattice2[2][1])
            layer_AorB_C[i][m][2] = (layer_A_or_B_D[i][m][0] * lattice1[0][2] + layer_A_or_B_D[i][m][1] \
                       * lattice1[1][2] + layer_A_or_B_D[i][m][2] * lattice2[2][2])
            layer_AorB_C[i][m][3] = (layer_A_or_B_D[i][m][3])
    return layer_AorB_C
def cal_dis1(x, y): 
    dist = ((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2) ** 0.5
    return dist
    

def retained_useful_atom(layer_A_or_B_C, bond_layer_A1_or_B1):  
    layer_A_or_B_C1 = copy.deepcopy(layer_A_or_B_C)
    for i in range(len(layer_A_or_B_C1) - 1):
        for m in layer_A_or_B_C1[i][:]:
            for n in layer_A_or_B_C1[i + 1][:]:
                bo = cal_dis1(m, n)
                if m in layer_A_or_B_C1[i]:
                    if n in layer_A_or_B_C1[i + 1]:
                        if (abs(bo - bond_layer_A1_or_B1[i]) < 0.0001): 
                            for l1 in layer_A_or_B_C1[i][:]:
                                if (m[3] == l1[3]):
                                    if (m[0] != l1[0] or m[1] != l1[1]):
                                        layer_A_or_B_C1[i].remove(l1)            
                            for l2 in layer_A_or_B_C1[i + 1][:]:
                                if n[3] == l2[3]:
                                    if (n[0] != l2[0] or n[1] != l2[1]):
                                        layer_A_or_B_C1[i + 1].remove(l2)
    layer_A_or_B_C2 = copy.deepcopy(layer_A_or_B_C1)
    for i in range(len(layer_A_or_B_C2)):
        for m in layer_A_or_B_C1[i][:]:
            for n in layer_A_or_B_C1[i][:]:
                if m in layer_A_or_B_C1[i] and n in layer_A_or_B_C1[i]:
                    if m[3] == n[3]:
                        if m[0] != n[0] or m[1] != n[1]:
                            layer_A_or_B_C2[i].remove(n)
    return layer_A_or_B_C2


def test1(layer_A_or_B_C_fina, C_A_or_B): 
    if len(layer_A_or_B_C_fina)>1:
        for i in range(len(layer_A_or_B_C_fina)):
            if abs(layer_A_or_B_C_fina[i][0][2] - 0) < 0.005:
                if C_A_or_B<0:
                    if layer_A_or_B_C_fina[i][0][2] - layer_A_or_B_C_fina[i + 1][0][2] > 3:
                        for m in range(len(layer_A_or_B_C_fina[i])):
                            layer_A_or_B_C_fina[i][m][2] = C_A_or_B
                elif C_A_or_B>0:
                    if layer_A_or_B_C_fina[i][0][2] - layer_A_or_B_C_fina[i - 1][0][2] > 3:
                        for m in range(len(layer_A_or_B_C_fina[i])):
                            layer_A_or_B_C_fina[i][m][2] = C_A_or_B
            elif abs(layer_A_or_B_C_fina[i][0][2] - C_A_or_B) < 0.005:
                if C_A_or_B<0:
                    if layer_A_or_B_C_fina[i][0][2] - layer_A_or_B_C_fina[i - 1][0][2] > 3:
                        for n in range(len(layer_A_or_B_C_fina[i])):
                            layer_A_or_B_C_fina[i][n][2] = 0
                elif C_A_or_B>0:
                    if layer_A_or_B_C_fina[i][0][2] - layer_A_or_B_C_fina[i - 1][0][2] > 3:
                        for p in range(len(layer_A_or_B_C_fina[i])):
                            layer_A_or_B_C_fina[i][p][2] = 0
    return layer_A_or_B_C_fina

if len(layer_A1) > 1:
    bond_layer_A = cal_bond_layer(layer_A1)
    if max(bond_layer_A) > 4.5: 

        layer_A_D = cartesian_to_direct(layer_A1, la_A, la_A)
        layer_A_D1 = copy.deepcopy(layer_A_D)
        for i in range(len(layer_A_D)): 
            for m in range(len(layer_A_D[i])):
                if abs(layer_A_D[i][m][0] - 0) < 0.0001:
                    layer_A_D1[i].append(list(layer_A_D[i][m]))
                    layer_A_D1[i][-1][0] = 1
                elif abs(layer_A_D[i][m][0] - 1) < 0.0001:
                    layer_A_D1[i].append(list(layer_A_D[i][m]))
                    layer_A_D1[i][-1][0] = 0
        for i in range(len(layer_A_D)):
            for m in range(len(layer_A_D[i])):
                if abs(layer_A_D[i][m][1] - 0) < 0.0001:
                    layer_A_D1[i].append(list(layer_A_D[i][m]))
                    layer_A_D1[i][-1][1] = 1
                elif abs(layer_A_D[i][m][1] - 1) < 0.0001:
                    layer_A_D1[i].append(list(layer_A_D[i][m]))
                    layer_A_D1[i][-1][1] = 0
        layer_A_C = direct_to_cartesian(layer_A_D1, la_A, la_A)
        bond_layer_A1 = cal_bond_layer(layer_A_C)
        layer_A_C_fina = retained_useful_atom(layer_A_C, bond_layer_A1) 
        layer1 = test1(layer_A_C_fina, C_A) 
    else:
        layer1 = test1(layer_A, C_A)
else:
    layer1 = test1(layer_A, C_A)
            
if len(layer_B1) > 1:
    bond_layer_B = cal_bond_layer(layer_B1)
    if max(bond_layer_B) > 4.5:
        layer_B_D = cartesian_to_direct(layer_B1, la_B, la_B)
        layer_B_D1 = copy.deepcopy(layer_B_D)
        for i in range(len(layer_B_D)): 
            for m in range(len(layer_B_D[i])):
                if abs(layer_B_D[i][m][0] - 0) < 0.0001:
                    layer_B_D1[i].append(list(layer_B_D[i][m]))
                    layer_B_D1[i][-1][0] = 1
                elif abs(layer_B_D[i][m][0] - 1) < 0.0001:
                    layer_B_D1[i].append(list(layer_B_D[i][m]))
                    layer_B_D1[i][-1][0] = 0
        for i in range(len(layer_B_D)):
            for m in range(len(layer_B_D[i])):
                if abs(layer_B_D[i][m][1] - 0) < 0.0001:
                    layer_B_D1[i].append(list(layer_B_D[i][m]))
                    layer_B_D1[i][-1][1] = 1
                elif abs(layer_B_D[i][m][1] - 1) < 0.0001:
                    layer_B_D1[i].append(list(layer_B_D[i][m]))
                    layer_B_D1[i][-1][1] = 0
        layer_B_C = direct_to_cartesian(layer_B_D1, la_B, la_B)
        bond_layer_B1 = cal_bond_layer(layer_B_C)
        layer_B_C_fina = retained_useful_atom(layer_B_C, bond_layer_B1)
        layer2 = test1(layer_B_C_fina, C_B)
    else:
        layer2 = test1(layer_B, C_B)
else:
    layer2 = test1(layer_B, C_B)


def test2(layerA_or_B, C_A_or_B): 
    inter_spac = []
    if len(layerA_or_B) > 1:
        for i in range(len(layerA_or_B) - 1):
            inter_spac.append(abs(layerA_or_B[i][0][2] - layerA_or_B[i + 1][0][2]))
    if C_A_or_B > 0:
        inter_spac.append(abs(layerA_or_B[-1][0][2] - 0) + abs(C_A_or_B - layerA_or_B[0][0][2]))
    if C_A_or_B < 0:
        inter_spac.append(abs(layerA_or_B[-1][0][2] - C_A_or_B) + abs(0 - layerA_or_B[0][0][2]))
    return inter_spac
        
def writen(layer1_or_2, atos_ele_C_B_or_A):
    for i in range(len(layer1_or_2)):
        for m in range(len(layer1_or_2[i])):
            for n in range(len(atos_ele_C_B_or_A)):
                if layer1_or_2[i][m][3] == atos_ele_C_B_or_A[n][3]:
                    atos_ele_C_B_or_A[n][2] = layer1_or_2[i][m][2]
                    atos_ele_C_B_or_A[n][1] = layer1_or_2[i][m][1]
                    atos_ele_C_B_or_A[n][0] = layer1_or_2[i][m][0]
    return atos_ele_C_B_or_A
atoms_A = writen(layer1, atos_ele_C_A)
atoms_B = writen(layer2, atos_ele_C_B)

POS1 = 'POSCAR111'
POS2 = 'POSCAR222'
fractor = 1.0

def save(fractor_A_or_B, POS1_or_2, la_A_or_B, ele_A_or_B, atom_num_A_or_B, atoms_A_or_B):
    title = "heter"
    scaling = str(fractor_A_or_B)
    fp4 = open(POS1_or_2, mode='w')
    float_width = 25
    str_width = 5
    fp4.write( title )
    fp4.write("\n")
    fp4.write( scaling )
    fp4.write("\n")
    fp4.write(str(la_A_or_B[0][0]).rjust(float_width, ' '))
    fp4.write(str(la_A_or_B[0][1]).rjust(float_width, ' '))
    fp4.write(str(la_A_or_B[0][2]).rjust(float_width, ' '))
    fp4.write("\n")
    fp4.write(str(la_A_or_B[1][0]).rjust(float_width, ' '))
    fp4.write(str(la_A_or_B[1][1]).rjust(float_width, ' '))
    fp4.write(str(la_A_or_B[1][2]).rjust(float_width, ' '))
    fp4.write("\n")
    fp4.write(str(la_A_or_B[2][0]).rjust(float_width, ' '))
    fp4.write(str(la_A_or_B[2][1]).rjust(float_width, ' '))
    fp4.write(str(la_A_or_B[2][2]).rjust(float_width, ' '))
    fp4.write("\n")
    for m in range(len(ele_A_or_B)):
        fp4.write(str(ele_A_or_B[m]).rjust(str_width, ' '))
    fp4.write("\n")
    for i in range(len(atom_num_A_or_B)):
        fp4.write(str(atom_num_A_or_B[i]).rjust(str_width, ' '))
    fp4.write("\n")
    fp4.write( "Cartesian" )
    fp4.write("\n")
    for i in range(len(atoms_A_or_B)):
        fp4.write(str(atoms_A_or_B[i][0]).rjust(float_width, ' '))
        fp4.write(str(atoms_A_or_B[i][1]).rjust(float_width, ' '))
        fp4.write(str(atoms_A_or_B[i][2]).rjust(float_width, ' '))
        fp4.write(str(atoms_A_or_B[i][3]).rjust(str_width, ' '))
        fp4.write("\n")
    fp4.close()
save(fractor, POS1, la_A, ele_A, atom_num_A, atoms_A)
save(fractor, POS2, la_B, ele_B, atom_num_B, atoms_B)

[fractor_A, la_A, ele_A, atom_num_A, atos_ele_C_A, C_A, atos_ele_C_sort_A] = read_first("POSCAR111")
[fractor_B, la_B, ele_B, atom_num_B, atos_ele_C_B, C_B, atos_ele_C_sort_B] = read_first("POSCAR222")
layer_A = dividing_layer(atos_ele_C_sort_A)
layer_B = dividing_layer(atos_ele_C_sort_B)
def cut(inter_spac_A_or_B, layer1_or_2): 
    if len(layer1_or_2) > 1: 
        for i in range(len(layer1_or_2) - 1):
            if inter_spac_A_or_B[i] > 3.5:
                break
        for m in range(len(layer1_or_2) - i - 1):
            del layer1_or_2[3]
        ele1 = []
        for i in range(len(layer1_or_2)):
            for m in range(len(layer1_or_2[i])):
                ele1.append(layer1_or_2[i][m][3][:-1])
        ele1_set = [[],[]]
        ele1_set[0] = set(ele1)
        for i in ele1_set[0]:
            ele1_set[1].append(ele1.count(i))
    return layer1_or_2, ele1_set

if len(layer1) > 1:
    inter_spac_A = test2(layer_A, C_A)
    if len([x for x in inter_spac_A if x > 3.5]) > 1: 
        print("Found that the first structure is the two-dimensional structure of the bulk")
        sele1 = 'y'
        if sele1 == 'y':
            print("Now start cutting surfaces...")
            [layer11, ele_set_A] = cut(inter_spac_A, layer1)
            ele_A = []
            for i in ele_set_A[0]:
                ele_A.append(i)
            atom_A = []
            for i in range(len(ele_A)):
                for m in range(len(layer11)):
                    for n in range(len(layer11[m])):
                        if layer11[m][n][3][:-1] == ele_A[i]:
                            atom_A.append(layer11[m][n])
            save(fractor_A, POS1, la_A, ele_A, ele_set_A[1], atom_A)
    elif len([x for x in inter_spac_A if x > 3.5]) == 1:
        for i in range(len(inter_spac_A) - 1):
            if inter_spac_A[i] > 3.5:
                for m in range(i + 1):
                    for n in range(len(layer_A[m])):
                        layer_A[m][n][2] = layer_A[m][n][2] - abs(C_A)
                for i in range(len(layer_A)):
                    for m in range(len(layer_A[i])):
                        layer_A[i][m][2] = layer_A[i][m][2] + abs(C_A) / 2
        at_A = writen(layer_A, atos_ele_C_A)
        save(fractor_A, POS1, la_A, ele_A, atom_num_A, at_A)
if len(layer2) > 1:
    inter_spac_B = test2(layer_B, C_B)
    if len([x for x in inter_spac_B if x > 3.5]) > 1: 
        print("Found that the second structure is the two-dimensional structure of the bulk")
        sele2 = 'y'
        if sele2 == 'y':
            print("Now start cutting surfaces...")
            [layer22, ele_set_B] = cut(inter_spac_B, layer2)
            ele_B = []
            for i in ele_set_B[0]:
                ele_B.append(i)
            atom_B = []
            for i in range(len(ele_B)):
                for m in range(len(layer22)):
                    for n in range(len(layer22[m])):
                        if layer22[m][n][3][:-1] == ele_B[i]:
                            atom_B.append(layer22[m][n])
            save(fractor_B, POS2, la_B, ele_B, ele_set_B[1], atom_B)
    elif len([x for x in inter_spac_B if x > 3.5]) == 1:
        for i in range(len(inter_spac_B) - 1):
            if inter_spac_B[i] > 3.5:
                for m in range(i + 1):
                    for n in range(len(layer_B[m])):
                        layer_B[m][n][2] = layer_B[m][n][2] - abs(C_B)
                for i in range(len(layer_B)):
                    for m in range(len(layer_B[i])):
                        layer_B[i][m][2] = layer_B[i][m][2] + abs(C_B) / 2
        at_B = writen(layer_B, atos_ele_C_B)
        save(fractor_B, POS2, la_B, ele_B, atom_num_B, at_B)

