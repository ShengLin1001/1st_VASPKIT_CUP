#! /bin/sh
echo "*************************************************************"
echo "*************************************************************"
echo "Author:JunLuo and JiZhou"
echo "Email:junluo@hnu.edu.cn  Hunan University"
echo "Function: A program that can be used to build heterostructure with two POSCAR"
echo "Version:1.0"
echo "First released time: 2019.07.30"
echo "Last modified time: 2019.07.30"
echo "*************************************************************"
echo "*************************************************************"
cp $1 POSCAR
echo 602 | vaspkit | grep 123
echo "The target file has been converted into primitive cell!"
mv PRIMCELL.vasp POSCAR1


cp $2 POSCAR
echo 602 | vaspkit | grep 123
echo "The target file has been converted into primitive cell!"
mv PRIMCELL.vasp POSCAR2

python3 find_hetero.py

