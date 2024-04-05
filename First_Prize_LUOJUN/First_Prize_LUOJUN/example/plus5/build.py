#!/usr/bin/env python3
from numpy import *
import numpy as np
import math
import os
from fractions import Fraction

os.system('python3 pe-building.py')

os.system('cp POSCAR111 POSCAR')
os.system('(echo 406;echo 1;echo 4062) | vaspkit | grep 123')
os.system('mv POSCAR_REV POSCAR11')

os.system('cp POSCAR222 POSCAR')
os.system('(echo 406;echo 1;echo 4062) | vaspkit | grep 123')
os.system('mv POSCAR_REV POSCAR22')

os.system('cp POSCAR11 POSCAR')
os.system('(echo 406;echo 1;echo 4061) | vaspkit | grep 123')
os.system('mv POSCAR_REV POSCAR11_C')

os.system('cp POSCAR22 POSCAR')
os.system('(echo 406;echo 1;echo 4061) | vaspkit | grep 123')
os.system('mv POSCAR_REV POSCAR22_C')

os.system('cp POSCAR11 POSCAR11_D')
os.system("cat POSCAR11_D | sed -n '9,$p'| awk '{print $1,$2,$3,$4}' > POSCAR11_D_atom")
os.system("cat POSCAR11_C | sed -n '9,$p'| awk '{print $1,$2,$3,$4}' > POSCAR11_C_atom")
os.system("cat POSCAR11_D | sed -n '3,5p'| awk '{print $1,$2,$3}' > POSCAR11_D_lattice")

os.system('cp POSCAR22 POSCAR22_D')
os.system("cat POSCAR22_D | sed -n '9,$p'| awk '{print $1,$2,$3,$4}' > POSCAR22_D_atom")
os.system("cat POSCAR22_C | sed -n '9,$p'| awk '{print $1,$2,$3,$4}' > POSCAR22_C_atom")
os.system("cat POSCAR22_D | sed -n '3,5p'| awk '{print $1,$2,$3}' > POSCAR22_D_lattice")

print("Successfully extract the necessary information.")


def read_lattice_and_atom(POS, lattice, atom, atom_C):
    fp1 = open(POS, mode='r')
    fp1.readline()
    fractor = float(fp1.readline())
    fp1.readline()
    fp1.readline()
    fp1.readline()
    ele = fp1.readline()
    ele = ele.strip('\n').split()
    atom_num = list(fp1.readline().split())
    atom_num = list(map(int, atom_num))
    fp1.close()
    
    fp2 = open(lattice, mode='r')
    lats = fp2.readlines()    
    A = zeros((3,3),dtype=float)
    A_row = 0
    for lat in lats:              
        list1 = lat.strip('\n').split(' ')      
        A[A_row:] = list1[0:3]                   
        A_row+=1          
    A = fractor * A
    angle = math.acos((A[0,0] * A[1,0] + A[0,1] * A[1,1]) / (
                ((A[0,0] ** 2 + A[0,1] ** 2) ** 0.5) * (
                    (A[1,0] ** 2 + A[1,1] ** 2) ** 0.5)))
    a = ((A[0,0] ** 2 + A[0,1] ** 2 + A[0,2] ** 2) ** 0.5) 
    b = ((A[1,0] ** 2 + A[1,1] ** 2 + A[1,2] ** 2) ** 0.5) 
    fp2.close()
    
    fp3 = open(atom, mode='r')
    atos = fp3.readlines()
    row = len(atos)
    fp3.close()
    
    fp4 = open(atom_C, mode='r')
    atos_ele_C_sort = [[] for m in range(row)]
    for i in range(row):
        atos_ele_C_sort[i][:] = list(fp4.readline().split())
        atos_ele_C_sort[i][0:3] = list(map(float, atos_ele_C_sort[i][0:3]))
    atos_ele_C_sort.sort(key=lambda x: x[2], reverse=True)
    fp4.close()
    
    fp6 = open(POS,mode='r')
    atos_D = [[] for m in range(row)]
    fp6.readline()
    fp6.readline()
    fp6.readline()
    fp6.readline()
    fp6.readline()
    fp6.readline()
    fp6.readline()
    fp6.readline()
    for i in range(row):
        atos_D[i][:] = list(fp6.readline().split())
        atos_D[i][0:3] = list(map(float, atos_D[i][0:3]))
    fp6.close()
    
    fp5 = open(atom_C, mode='r')
    atos_ele_C = [[] for m in range(row)]
    for i in range(row):
        atos_ele_C[i][:] = list(fp5.readline().split())
        atos_ele_C[i][0:3] = list(map(float, atos_ele_C[i][0:3]))
    return angle, A, atos_ele_C, ele, atom_num, a, b, atos_ele_C_sort, atos_D
    fp5.close()

[angle1, lattice_A, atom_ele_C_A, ele_A, atom_num_A, lat_conA_a, lat_conA_b, atos_ele_C_sort_A, atos_D_A] = read_lattice_and_atom("POSCAR11", "POSCAR11_D_lattice", "POSCAR11_D_atom", "POSCAR11_C_atom")
[angle2, lattice_B, atom_ele_C_B, ele_B, atom_num_B, lat_conB_a, lat_conB_b, atos_ele_C_sort_B, atos_D_B] = read_lattice_and_atom("POSCAR22", "POSCAR22_D_lattice", "POSCAR22_D_atom", "POSCAR22_C_atom")
lat_con_aver_a = (lat_conA_a + lat_conB_a) / 2
lat_con_aver_b = (lat_conA_b + lat_conB_b) / 2
angle_aver = (angle1 + angle2) / 2
vector = [[] for i in range(3)]
vector[0] = [math.cos(angle_aver / 2) * lat_con_aver_a, math.sin(angle_aver / 2) * lat_con_aver_a, 0]
vector[1] = [math.cos(angle_aver / 2) * lat_con_aver_b, math.sin(angle_aver / 2) * lat_con_aver_b * (-1), 0]
vector[2] = [0, 0, 0]
#vector = np.mat(vector)
 
print("Now we are starting to build a heterostructure")
print("Please customize the lattice parameters of the heterostructure:")
print("1) Lattice parameters of the lower material")
print("2) Lattice parameters of the upper material")
print("3) Average of the lattice parameters of the upper and lower materials")
lattice_select = int(input("Please customize the lattice parameters of the heterostructure:"))
vac_thickness = float(input("Please customize the thickness of vacuum layer of the heterostructure (unit: Angstrom):"))
inter_spacing = float(input("Please customize the layer spacing of the heterostructure (unit: Angstrom):"))

def cal_dis(coor1, coor2, q, w): 
    distance = ((coor1[q][0] - coor2[w][0]) ** 2 + (coor1[q][1] - coor2[w][1]) ** 2 + (coor1[q][2] - coor2[w][2]) ** 2) ** 0.5
    return distance

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
layer_B = dividing_layer(atos_ele_C_sort_B) 
layer_A = dividing_layer(atos_ele_C_sort_A)

def cal_bond_length_1(atom_ele_C_A_or_B, atom_num_A_B0): 
    atom_dis1 = []
    for i in range(atom_num_A_B0 - 1):
        atom_dis1.append(cal_dis(atom_ele_C_A_or_B, atom_ele_C_A_or_B, 0, i + 1))
        bond_length1 = min(atom_dis1)
    return bond_length1

def cal_bond_length_2(layer_A_or_B):  
    atom_dis = []
    for m in range(len(layer_A_or_B[0])):
        for n in range(len(layer_A_or_B[1])):
            atom_dis.append(cal_dis(layer_A_or_B[0], layer_A_or_B[1], m, n))
    bond_length = min(atom_dis)
    return bond_length

if len(layer_A) == 1:
    bond_length_A = cal_bond_length_1(atom_ele_C_A, atom_num_A[0])
else:
    bond_length_A = cal_bond_length_2(layer_A)

if len(layer_B) == 1:
    bond_length_B = cal_bond_length_1(atom_ele_C_B, atom_num_B[0])
else:
    bond_length_B = cal_bond_length_2(layer_B)


def cal_dividing_layer_space(layer_A_or_B, bond_length_A_or_B):
    std_atom_infor = [0 for i in range(len(layer_A_or_B) - 1)]
    layer_std = [0 for i in range(len(layer_A_or_B) - 1)]
    for i in range(len(layer_A_or_B) - 1): 
        for m in range(len(layer_A_or_B[i])):
            for n in range(len(layer_A_or_B[i + 1])):
                space = cal_dis(layer_A_or_B[i], layer_A_or_B[i+1], m, n)
                if abs(space - bond_length_A_or_B) < 0.001:
                    layer_std[i] = layer_A_or_B[i][m]
                    std_atom_infor[i] = layer_A_or_B[i + 1][n]
    return layer_std, std_atom_infor
[layer_std_A, std_atom_infor_A] = cal_dividing_layer_space(layer_A, bond_length_A)
[layer_std_B, std_atom_infor_B] = cal_dividing_layer_space(layer_B, bond_length_B)
def direct_to_cartesian(atos_D_A_or_B, lattice1,lattice2): 
    atom1_latt2 = [[] for m in range(len(atos_D_A_or_B))]
    for i in range(len(atos_D_A_or_B)):
        atom1_latt2[i].append(atos_D_A_or_B[i][0] * lattice1[0][0] + atos_D_A_or_B[i][1] * lattice1[1][0] + atos_D_A_or_B[i][2] * lattice2[2][0])
        atom1_latt2[i].append(atos_D_A_or_B[i][0] * lattice1[0][1] + atos_D_A_or_B[i][1] * lattice1[1][1] + atos_D_A_or_B[i][2] * lattice2[2][1])
        atom1_latt2[i].append(atos_D_A_or_B[i][0] * lattice1[0][2] + atos_D_A_or_B[i][1] * lattice1[1][2] + atos_D_A_or_B[i][2] * lattice2[2][2])
        atom1_latt2[i].append(atos_D_A_or_B[i][3])
    return atom1_latt2
atomB_lattAB = direct_to_cartesian(atos_D_B, lattice_A,lattice_B) 
atomA_lattBA = direct_to_cartesian(atos_D_A, lattice_B,lattice_A) 

atomB_latt_he_B = direct_to_cartesian(atos_D_B, vector,lattice_B) 
atomA_latt_he_A = direct_to_cartesian(atos_D_A, vector,lattice_A) 

def cartesian_to_direct(atomB_or_A_he_C, lattice1, lattice2):
    transit = [[] for i in range(len(atomB_or_A_he_C))]
    for i in range(len(atomB_or_A_he_C)):
        transit[i].append(atomB_or_A_he_C[i][0]) 
        transit[i].append(atomB_or_A_he_C[i][1])
        transit[i].append(atomB_or_A_he_C[i][2]) 
#    atom_C = np.transpose(transit)
    lattice = [[] for i in range(len(lattice1))]
    for i in range(len(lattice1) - 1):
        lattice[i].append(lattice1[i][0])
        lattice[i].append(lattice1[i][1])
        lattice[i].append(lattice1[i][2])
    lattice[2].append(lattice2[2][0])
    lattice[2].append(lattice2[2][1])
    lattice[2].append(lattice2[2][2])
    lattice_tran_mat = np.mat(np.transpose(lattice))
    convmat = np.linalg.inv(lattice_tran_mat)
    latt_inv = convmat.tolist()
    atomB_or_A_he_D = [[] for m in range(len(atomB_or_A_he_C))]
    for i in range(len(transit)):
        atomB_or_A_he_D[i].append(latt_inv[0][0] * transit[i][0] + latt_inv[0][1] * transit[i][1] + latt_inv[0][2] * transit[i][2])
        atomB_or_A_he_D[i].append(latt_inv[1][0] * transit[i][0] + latt_inv[1][1] * transit[i][1] + latt_inv[1][2] * transit[i][2])
        atomB_or_A_he_D[i].append(latt_inv[2][0] * transit[i][0] + latt_inv[2][1] * transit[i][1] + latt_inv[2][2] * transit[i][2])
        atomB_or_A_he_D[i].append(atomB_or_A_he_C[i][3])
    return atomB_or_A_he_D

def tran_of_each_layer(layer_std_A_or_B, atom_lattAB, std_atom_infor_A_or_B, bond_length_A_or_B):
    for i in range(len(layer_std_A_or_B)):
        for m in range(len(atom_lattAB)):
            for n in range(len(atom_lattAB)):
                if layer_std_A_or_B[i][3] == atom_lattAB[m][3]:
                    if std_atom_infor_A_or_B[i][3] == atom_lattAB[n][3]:
                        atom_tran = atom_lattAB[m][2] - (bond_length_A_or_B **  \
                        2 - (atom_lattAB[m][0] - atom_lattAB[n][0]) ** 2 - (atom_lattAB[m][1] -  \
                        atom_lattAB[n][1]) ** 2) ** 0.5
                        tt = atom_lattAB[n][2]
                        for o in range(len(atom_lattAB)):
                            if atom_lattAB[o][2] == tt:
                                atom_lattAB[o][2] = atom_tran
    return atom_lattAB

atomB_he_C = tran_of_each_layer(layer_std_B, atomB_lattAB, std_atom_infor_B, bond_length_B)
atomB_he_D = cartesian_to_direct(atomB_he_C, lattice_A, lattice_B)
atomA_he_C = tran_of_each_layer(layer_std_A, atomA_lattBA, std_atom_infor_A, bond_length_A)
atomA_he_D = cartesian_to_direct(atomA_he_C, lattice_B, lattice_A)

atomB_hehe_C = tran_of_each_layer(layer_std_B, atomB_latt_he_B, std_atom_infor_B, bond_length_B)
atomB_hehe_D = cartesian_to_direct(atomB_hehe_C, vector, lattice_B)
atomA_hehe_C = tran_of_each_layer(layer_std_A, atomA_latt_he_A, std_atom_infor_A, bond_length_A)
atomA_hehe_D = cartesian_to_direct(atomA_hehe_C, vector, lattice_A)

for i in range(len(atos_D_A)):
    atos_D_A[i].pop(3)
for i in range(len(atos_D_B)):
    atos_D_B[i].pop(3)
for i in range(len(atomA_he_D)):
    atomA_he_D[i].pop(3)
for i in range(len(atomB_he_D)):
    atomB_he_D[i].pop(3)
for i in range(len(atomA_hehe_D)):
    atomA_hehe_D[i].pop(3)
for i in range(len(atomB_hehe_D)):
    atomB_hehe_D[i].pop(3)
    
if lattice_select == 1:
    atom_A = np.mat(atos_D_A)
    atom_B = np.mat(atomB_he_D)
elif lattice_select == 2:
    atom_A = np.mat(atomA_he_D)
    atom_B = np.mat(atos_D_B)
else:
    atom_A = np.mat(atomA_hehe_D)
    atom_B = np.mat(atomB_hehe_D)


dis_A = (atom_A[:,2].max(axis=0) - atom_A[:,2].min(axis=0)) * abs(lattice_A[2,2])  
dis_B = (atom_B[:,2].max(axis=0) - atom_B[:,2].min(axis=0)) * abs(lattice_B[2,2])
C = dis_A + dis_B + vac_thickness + inter_spacing
C_A = abs(lattice_A[2,2])
C_B = abs(lattice_B[2,2])
vector = np.mat(vector)
vector[2,2] = C
lattice_aver = vector
atom_A_new = atom_A
atom_B_new = atom_B
slab_pos = vac_thickness / 2
distance_A_Direct = atom_A[:,2].min(axis=0) - (slab_pos / C_A)
atom_A_new[:,2] = (atom_A[:,2] - distance_A_Direct) * C_A / C  
A1_max = atom_A_new[:,2].max(axis=0)   
B_min_Cartesian = (atom_B[:,2].min(axis=0)) * C_B
B1_min_Cartesian = A1_max * C + inter_spacing    
distance_B_Cartesian = B_min_Cartesian - B1_min_Cartesian
atom_B_new[:,2] = (atom_B[:,2] * C_B - distance_B_Cartesian) / C
lattice_A[2,2] = C
lattice_B[2,2] = C
name0 = "hetero"
name1 = str(lattice_select)
name2 = str(vac_thickness)
name3 = str(inter_spacing)
name4 = ".vasp"
name_t = name0 + '_' + name1 + '_' + name2 + '_' + name3 + name4
def building_heterostructure(lattice_he): 
    title = "Heterostructure generated by building_heterostructure"
    scaling = "1.0"
    fp4 = open(name_t, mode='w')
    float_width = 25
    str_width = 5
    fp4.write( title )
    fp4.write("\n")
    fp4.write( scaling )
    fp4.write("\n")
    fp4.write(str(lattice_he[0,0]).rjust(float_width, ' '))
    fp4.write(str(lattice_he[0,1]).rjust(float_width, ' '))
    fp4.write(str(lattice_he[0,2]).rjust(float_width, ' '))
    fp4.write("\n")
    fp4.write(str(lattice_he[1,0]).rjust(float_width, ' '))
    fp4.write(str(lattice_he[1,1]).rjust(float_width, ' '))
    fp4.write(str(lattice_he[1,2]).rjust(float_width, ' '))
    fp4.write("\n")
    fp4.write(str(lattice_he[2,0]).rjust(float_width, ' '))
    fp4.write(str(lattice_he[2,1]).rjust(float_width, ' '))
    fp4.write(str(lattice_he[2,2]).rjust(float_width, ' '))
    fp4.write("\n")
    for m in range(len(ele_A)):
        fp4.write(str(ele_A[m]).rjust(str_width, ' '))
    for n in range(len(ele_B)):
        fp4.write(str(ele_B[n]).rjust(str_width, ' '))
    fp4.write("\n")
    for i in range(len(atom_num_A)):
        fp4.write(str(atom_num_A[i]).rjust(str_width, ' '))
    for i in range(len(atom_num_B)):
        fp4.write(str(atom_num_B[i]).rjust(str_width, ' '))
    fp4.write("\n")
    fp4.write( "Direct" )
    fp4.write("\n")
    atom_A_new1 = atom_A_new.tolist()
    atom_B_new1 = atom_B_new.tolist()
    for i in range(len(atom_A_new1)):
        fp4.write(str(atom_A_new1[i][0]).rjust(float_width, ' '))
        fp4.write(str(atom_A_new1[i][1]).rjust(float_width, ' '))
        fp4.write(str(atom_A_new1[i][2]).rjust(float_width, ' '))
        fp4.write("\n")
    for i in range(len(atom_B_new1)):
        fp4.write(str(atom_B_new1[i][0]).rjust(float_width, ' '))
        fp4.write(str(atom_B_new1[i][1]).rjust(float_width, ' '))
        fp4.write(str(atom_B_new1[i][2]).rjust(float_width, ' '))
        fp4.write("\n")
    fp4.close()
if lattice_select == 1:
    building_heterostructure(lattice_A)
    angle_hehe = angle1
elif lattice_select == 2:
    building_heterostructure(lattice_B)
    angle_hehe = angle2
else:
    building_heterostructure(lattice_aver)
    angle_hehe = angle_aver
print("Written %s File ..." %(name_t))
print("Successfully build a heterostructure you want!")


print("Next you can also customize the horizontal offset of the upper or lower material.")
next_hor = input("Do you want to continue? (y/n): ")
if next_hor == 'y':
    print("Which layer of material do you want to translate?")
    print("1) The lower material")
    print("2) The upper material")
    tran_A_or_B = float(input())
    
    def hori_displace(vec_x, vec_y, angle_he):
        x_dis = vec_x / math.sin(angle_he) - vec_y / math.tan(angle_he)
        y_dis = -vec_x / math.tan(angle_he) + vec_y / math.sin(angle_he)
        return x_dis, y_dis

    print("Enter the translation vector :")
    print("(MUST be two fractions, e.g., 1/2 1/3)")
    aa1,bb1 = input().split()
    aa2 = str(round(float(Fraction(aa1)), 2))
    bb2 = str(round(float(Fraction(bb1)), 2))
    aa = -Fraction(aa1)
    bb = -Fraction(bb1)

    if abs(aa) > 1 or abs(bb) > 1:
        print("You must enter two numbers whose absolute values are less than 1!")
    else:
        [xx_dis, yy_dis] = hori_displace(aa, bb, angle_hehe)
        
        name00 = "heter_dised"
        if tran_A_or_B == 1:
            name01 = "lower"
        elif tran_A_or_B == 2:
            name01 = "upper"
        name_t_dised = name00 + '_' + name1 + '_' + name2 + '_' + name3 + '_' + name01 +  \
                        '_' + aa2 + '_' + bb2 + name4
        def building_hetero_hor(at_dis, lattice_he):
            atom_dised = zeros((len(at_dis),3),dtype=float)
            for i in range(len(at_dis)):
                atom_dised[i,0] =  at_dis[i,0] + xx_dis
                atom_dised[i,1] =  at_dis[i,1] + yy_dis
                atom_dised[i,2] =  at_dis[i,2]
            for i in range(len(atom_dised)):
                if atom_dised[i,0] > 1:
                    atom_dised[i,0] = atom_dised[i,0] - 1
                elif atom_dised[i,0] < -1:
                    atom_dised[i,0] = atom_dised[i,0] + 1
            for i in range(len(atom_dised)):
                if atom_dised[i,1] > 1:
                    atom_dised[i,1] = atom_dised[i,1] - 1
                elif atom_dised[i,1] < -1:
                    atom_dised[i,1] = atom_dised[i,1] + 1
            title = "Heter_dis generated by building_hetero_hor"
            scaling = "1.0"
            fp4 = open(name_t_dised, mode='w')
            float_width = 25
            str_width = 5
            fp4.write( title )
            fp4.write("\n")
            fp4.write( scaling )
            fp4.write("\n")
            fp4.write(str(lattice_he[0,0]).rjust(float_width, ' '))
            fp4.write(str(lattice_he[0,1]).rjust(float_width, ' '))
            fp4.write(str(lattice_he[0,2]).rjust(float_width, ' '))
            fp4.write("\n")
            fp4.write(str(lattice_he[1,0]).rjust(float_width, ' '))
            fp4.write(str(lattice_he[1,1]).rjust(float_width, ' '))
            fp4.write(str(lattice_he[1,2]).rjust(float_width, ' '))
            fp4.write("\n")
            fp4.write(str(lattice_he[2,0]).rjust(float_width, ' '))
            fp4.write(str(lattice_he[2,1]).rjust(float_width, ' '))
            fp4.write(str(lattice_he[2,2]).rjust(float_width, ' '))
            fp4.write("\n")
            for m in range(len(ele_A)):
                fp4.write(str(ele_A[m]).rjust(str_width, ' '))
            for n in range(len(ele_B)):
                fp4.write(str(ele_B[n]).rjust(str_width, ' '))
            fp4.write("\n")
            for i in range(len(atom_num_A)):
                fp4.write(str(atom_num_A[i]).rjust(str_width, ' '))
            for i in range(len(atom_num_B)):
                fp4.write(str(atom_num_B[i]).rjust(str_width, ' '))
            fp4.write("\n")
            fp4.write( "Direct" )
            fp4.write("\n")
            if tran_A_or_B == 1:
                atom_A_new1 = atom_dised.tolist()
                atom_B_new1 = atom_B_new.tolist()
            else:
                atom_A_new1 = atom_A_new.tolist()
                atom_B_new1 = atom_dised.tolist()
            for i in range(len(atom_A_new1)):
                fp4.write(str(atom_A_new1[i][0]).rjust(float_width, ' '))
                fp4.write(str(atom_A_new1[i][1]).rjust(float_width, ' '))
                fp4.write(str(atom_A_new1[i][2]).rjust(float_width, ' '))
                fp4.write("\n")
            for i in range(len(atom_B_new1)):
                fp4.write(str(atom_B_new1[i][0]).rjust(float_width, ' '))
                fp4.write(str(atom_B_new1[i][1]).rjust(float_width, ' '))
                fp4.write(str(atom_B_new1[i][2]).rjust(float_width, ' '))
                fp4.write("\n")
            fp4.close()
        
        if tran_A_or_B == 1:
            if lattice_select == 1:
                building_hetero_hor(atom_A_new, lattice_A)
            elif lattice_select == 2:
                building_hetero_hor(atom_A_new, lattice_B)
            else:
                building_hetero_hor(atom_A_new, lattice_aver)
        else:
            if lattice_select == 1:
                building_hetero_hor(atom_B_new, lattice_A)
            elif lattice_select == 2:
                building_hetero_hor(atom_B_new, lattice_B)
            else:
                building_hetero_hor(atom_B_new, lattice_aver)
    print("Written %s File ..." %(name_t_dised))
    print("Successfully build another stacking style heterostructure you want!")

os.system('rm POSCAR11_* POSCAR22_*')








    


