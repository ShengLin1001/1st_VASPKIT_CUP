# -*- coding:utf-8 -*-

"""
@author: Nan-nan Zhang
@file: heteroJ
@time: 2019/07/22 19:39
"""

from pymatgen import Structure
from pymatgen.io.vasp.outputs import Poscar
from pymatgen.io.cif import CifWriter
import numpy as np


mismatchA = 0.001
mismatchL = 0.05
file1 = "04_Ni111_POSCAR"
file2 = "04_MoS2_POSCAR"
vacuum = 15
distance = 3

print("\n########### Begin ############")
print("Author: Nan-nan Zhang\n",
      "zhangnn18@mails.tsinghua.edu.cn\n",
      "First released time: 2019/07/23\n",
      "Last modified time: 2019/07/23\n",
      "Version: v0.01")
print("Functions: \nDo Heterostructures Junction build from two POSCAR\n",
      "Usage:\n",
      "Set the mismatch Angel, mismatch Length at first lines of scripts\n",
      "Input the file name file1 and file2\n",
      "Set the Vacuum and distance at first lines of scripts.")
print("##############################\n")


def readstruct(file,supercell=[1,1,1]):

    poscar1 = Poscar.from_file(file)    # 调用pymatgen读取POSCAR
    structure = poscar1.structure                   # 整理成pymatgen.core.structure类的变量
    structure.make_supercell(supercell)
    #print("output lattice parameters:")
    #print(structure.lattice)                        # 打印lattice参数(3*3 matrix)
    strucPara = str(structure.get_sorted_structure()).split('\n')    # 提取夹角和边长

    #print(file, "\t parameters:")
    #print(strucPara[2])
    #print(strucPara[3])
    return strucPara[2], strucPara[3]


def makeJunction(file1, file2, supercell1, supercell2, vacuum=15, distance=2):
    poscar1 = Poscar.from_file(file1)    # 调用pymatgen读取POSCAR
    structure1 = poscar1.structure                   # 整理成pymatgen.core.structure类的变量
    structure1.make_supercell(supercell1)
    poscar2 = Poscar.from_file(file2)    # 调用pymatgen读取POSCAR
    structure2 = poscar2.structure                   # 整理成pymatgen.core.structure类的变量
    structure2.make_supercell(supercell2)

    # find the biggest z, find the smallest z
    zlist1 = []
    for i in range(len(structure1.species)):
        zlist1.append(float(structure1.frac_coords[i][-1]))
    zmax1,zmin1 = max(zlist1), min(zlist1)
    thick1 = zmax1-zmin1
    zlist2 = []
    for i in range(len(structure2.species)):
        zlist2.append(float(structure2.frac_coords[i][-1]))
    zmax2,zmin2 = max(zlist2), min(zlist2)
    thick2 = zmax2 - zmin2


    abc1, angle1 = readstruct(file1, supercell1)
    abc1_clean = make_clean(abc1)

    thicktot = (thick1 + thick2)*float(abc1_clean[-1]) + distance + vacuum

    for i in range(len(structure2.species)):
        newcoordlist = []
        #print(structure2.frac_coords[i])

        newcoordlist.append(structure2.frac_coords[i][0])
        newcoordlist.append(structure2.frac_coords[i][1])
        newcoordlist.append(structure2.frac_coords[i][-1] + zmax1 - zmin2 + distance/float(abc1_clean[-1]))
        #print(structure2.frac_coords[i],newcoordlist)
        structure1.append(species=structure2.species[i], coords=newcoordlist, coords_are_cartesian=False)
    print("\n##############################")
    print("output " + "POSCAR" + file1 + file2)
    print("##############################\n")
    open("POSCAR" + file1 + file2, "w").write(str(Poscar(structure1).get_string(direct=False)))
    # modify z to vacuum thick
    f = open("POSCAR" + file1 + file2,'r',encoding='utf-8')
    f_new = open("POSCAR" + file1 + file2 + "withVaccume", 'w', encoding='utf-8')
    count = 0
    for line in f:
        count += 1
        if count == 5:
            line = "0.0 0.0 " + str(thicktot) + "\n"
            f_new.write(line)
            continue
        f_new.write(line)
    f.close()
    f_new.close()


def make_clean(abc):
    abc_clean = abc.strip().split(" ")
    while '' in abc_clean:
        abc_clean.remove('')
    return list(map(lambda x: float(x), abc_clean[-3:]))


supercell1 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
supercell2 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
print("\nread structure 1")
abc1, angle1 = readstruct(file1,supercell1)
print("\nread structure 2")
abc2, angle2 = readstruct(file2,supercell2)

abc1_clean = make_clean(abc1)
angle1_clean = make_clean(angle1)
abc2_clean = make_clean(abc2)
angle2_clean = make_clean(angle2)

print("\n##############################")
print("Do angle check")
mismatch_angle3 = abs(angle1_clean[2]-angle2_clean[2])/min((angle1_clean[2],angle2_clean[2]))
print("mismatch angle3:", mismatch_angle3)
print("##############################\n")

if mismatch_angle3 < mismatchA:
    print('Good! mismatch_angle statisfied.')
else:
    print('Oh! mismatch_angle NOT statisfied. Finding Supercell. Please wait a few second')
    judgeSmallest = 100
    for a in range(-2, 3):
        for b in range(-2, 3):
            for i in range(-2, 3):
                for j in range(-2, 3):
                    if a == 0 or i == 0:
                        continue

                    abc1, angle1 = readstruct(file1, [[a, b, 0],[0, 1, 0],[0, 0, 1]])
                    abc2, angle2 = readstruct(file2, [[i, j, 0],[0, 1, 0],[0, 0, 1]])
                    abc1_clean = make_clean(abc1)
                    angle1_clean = make_clean(angle1)
                    abc2_clean = make_clean(abc2)
                    angle2_clean = make_clean(angle2)
                    mismatch_angle3 = abs(angle1_clean[2] - angle2_clean[2]) / min(
                        (angle1_clean[2], angle2_clean[2]))
                    if mismatch_angle3 < mismatchA:
                        if abs(a)+abs(b)+abs(i)+abs(j) >= judgeSmallest:
                            continue
                        judgeSmallest = abs(a)+abs(b)+abs(i)+abs(j)
                        supercell1 = [[a, b, 0],[0, 1, 0],[0, 0, 1]]
                        supercell2 = [[i, j, 0],[0, 1, 0],[0, 0, 1]]
                        print("\nmismatch =", mismatch_angle3)
                        print("Find a new supercellgroup:", supercell1, supercell2)
                        print(abc1, angle1)
                        print(abc2, angle2)
                        break
                    else:
                        continue

print("\n##############################")
print("Angle match Use",supercell1,supercell2)
print("##############################\n")
abc1, angle1 = readstruct(file1, supercell1)
abc2, angle2 = readstruct(file2, supercell2)
abc1_clean = make_clean(abc1)
angle1_clean = make_clean(angle1)
abc2_clean = make_clean(abc2)
angle2_clean = make_clean(angle2)


mismatch_abc1 = abs(abc1_clean[0]-abc2_clean[0])/min(abs(abc2_clean[0]),abs(abc2_clean[0]))
mismatch_abc2 = abs(abc1_clean[1]-abc2_clean[1])/min(abs(abc2_clean[1]),abs(abc2_clean[1]))
print("\n##############################")
print("Do Length check")
print("mismatch abc1:", mismatch_abc1)
print("mismatch abc2:", mismatch_abc2)
print("##############################\n")

if mismatch_abc1 < mismatchL and mismatch_abc2 < mismatchL:
    print('Good! mismatch_Length statisfied.')
else:
    print('Oh! mismatch_Length NOT statisfied. Finding Supercell. Please wait a few second')
    judgeSmallest = 100
    for a in range(1, 9):
        for b in range(1, 9):
            for i in range(1, 9):
                for j in range(1, 9):
                    mismatch_abc1 = abs(abc1_clean[0]*a-abc2_clean[0]*i)/min(abs(abc2_clean[0])*a,abs(abc2_clean[0]*i))
                    mismatch_abc2 = abs(abc1_clean[1]*b-abc2_clean[1]*j)/min(abs(abc2_clean[1])*b,abs(abc2_clean[1]*j))
                    if mismatch_abc1 < mismatchL and mismatch_abc2 < mismatchL:
                        if abs(a) + abs(b) + abs(i) + abs(j) >= judgeSmallest:
                            continue
                        judgeSmallest = abs(a) + abs(b) + abs(i) + abs(j)
                        print("\nFind a supercellgroup:")
                        print("mismatch abc1:", mismatch_abc1)
                        print("mismatch abc2:", mismatch_abc2)
                        Lengthsuperlist = [a,b,i,j]
    if 'Lengthsuperlist' not in dir():
        print("@@@@@@@@@ reverse one of the vector and retry @@@@@@@@@@")
        judgeSmallest = 100
        for a in range(1, 9):
            for b in range(1, 9):
                for i in range(1, 9):
                    for j in range(1, 9):
                        mismatch_abc1 = abs(abc1_clean[1] * a - abc2_clean[0] * i) / min(abs(abc2_clean[1]) * a,
                                                                                         abs(abc2_clean[0] * i))
                        mismatch_abc2 = abs(abc1_clean[0] * b - abc2_clean[1] * j) / min(abs(abc2_clean[0]) * b,
                                                                                         abs(abc2_clean[1] * j))
                        if mismatch_abc1 < mismatchL and mismatch_abc2 < mismatchL:
                            if abs(a) + abs(b) + abs(i) + abs(j) >= judgeSmallest:
                                continue
                            judgeSmallest = abs(a) + abs(b) + abs(i) + abs(j)
                            print("\nFind a supercellgroup:")
                            print("mismatch abc1:", mismatch_abc1)
                            print("mismatch abc2:", mismatch_abc2)
                            Lengthsuperlist = [b, a, i, j]
                            print(Lengthsuperlist)

    supercell1[0][0] = supercell1[0][0] * Lengthsuperlist[0]
    supercell1[0][1] = supercell1[0][1] * Lengthsuperlist[0]
    supercell1[1][0] = supercell1[1][0] * Lengthsuperlist[1]
    supercell1[1][1] = supercell1[1][1] * Lengthsuperlist[1]

    supercell2[0][0] = supercell2[0][0] * Lengthsuperlist[2]
    supercell2[0][1] = supercell2[0][1] * Lengthsuperlist[2]
    supercell2[1][0] = supercell2[1][0] * Lengthsuperlist[3]
    supercell2[1][1] = supercell2[1][1] * Lengthsuperlist[3]

print("\n##############################")
print("Angel and Length match Use",supercell1,supercell2)
print("##############################\n")

print("\nread new structure 1")
abc1, angle1 = readstruct(file1,supercell1)
print(abc1, "\n", angle1)
print("\nread new structure 2")
abc2, angle2 = readstruct(file2,supercell2)
print(abc2, "\n", angle2)
abc1_clean = make_clean(abc1)
angle1_clean = make_clean(angle1)
abc2_clean = make_clean(abc2)
angle2_clean = make_clean(angle2)

makeJunction(file1, file2, supercell1, supercell2, vacuum, distance)

print("############## Summary ################")
print("The final supercell for ", file1)
print(supercell1[0],"\n",supercell1[1],"\n",supercell1[2],"\n")
print("The final supercell for ", file2)
print(supercell2[0],"\n",supercell2[1],"\n",supercell2[2],"\n")
print("Final mismatch of Angel:", abs(angle1_clean[2]-angle2_clean[2])/min((angle1_clean[2],angle2_clean[2])))
mismatch_abc1 = abs(abc1_clean[0]-abc2_clean[0])/min(abs(abc2_clean[0]),abs(abc2_clean[0]))
mismatch_abc2 = abs(abc1_clean[1]-abc2_clean[1])/min(abs(abc2_clean[1]),abs(abc2_clean[1]))
print("Final mismatch of Length a and b:", mismatch_abc1, " ", mismatch_abc2)
print("Use ", file1, " as substrate\n"
      "Distance between two slab is", distance, "Angstrom"
      "\nVacuum between two slab is", vacuum, "Angstrom")
print("############## Summary ################")