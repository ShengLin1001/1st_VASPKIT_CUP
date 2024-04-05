#!/bin/sh

echo "1) sorted_by_atomic_number.txt"
echo "2) sorted_by_custom_angle.txt"
echo "------------>>"
read -p "Which file data do you want to use:" file_number

if [ "$file_number" = 1 ]; then
    echo 'reading sorted_by_atomic_number.txt ....'
    read -p "Which row of data do you want to use in sorted_by_atomic_number.txt:" row
    x1=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $1}')
    y1=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $2}')
    x2=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $3}')
    y2=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $4}')
    z0=$(expr $x1 \* $y2 - $x2 \* $y1) 
    if [ "$z0" -gt 0 ]
    then
        z1=1
    else
        z1=-1
    fi
    xx1=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $5}')
    yy1=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $6}')
    xx2=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $7}')
    yy2=$(sed -n ${row}p sorted_by_atomic_number.txt | awk '{print $8}')
    zz0=$(expr $xx1 \* $yy2 - $xx2 \* $yy1)
    if [ "$zz0" -gt 0 ]
    then
        z2=1
    else
        z2=-1
    fi
    cp POSCAR1 POSCAR
    (echo 400; echo $x1 $y1 0; echo $x2 $y2 0; echo 0 0 $z1) | vaspkit | grep 123
    mv SUPERCELL.vasp POSCAR11
    cp POSCAR2 POSCAR
    (echo 400; echo $xx1 $yy1 0; echo $xx2 $yy2 0; echo 0 0 $z2) | vaspkit | grep 123
    mv SUPERCELL.vasp POSCAR22

else
    echo 'reading sorted_by_custom_angle.txt ....'
    read -p "Which row of data do you want to use in sorted_by_custom_angle.txt:" row
    x1=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $1}')
    y1=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $2}')
    x2=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $3}')
    y2=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $4}')
    z0=$(expr $x1 \* $y2 - $x2 \* $y1)
    if [ "$z0" -gt 0 ]
    then
        z1=1
    else
        z1=-1
    fi
    xx1=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $5}')
    yy1=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $6}')
    xx2=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $7}')
    yy2=$(sed -n ${row}p sorted_by_custom_angle.txt | awk '{print $8}')
    zz0=$(expr $xx1 \* $yy2 - $xx2 \* $yy1)
    if [ "$zz0" -gt 0 ]
    then
        z2=1
    else
        z2=-1
    fi
    cp POSCAR1 POSCAR
    (echo 400; echo $x1 $y1 0; echo $x2 $y2 0; echo 0 0 $z1) | vaspkit | grep 123
    mv SUPERCELL.vasp POSCAR11
    cp POSCAR2 POSCAR
    (echo 400; echo $xx1 $yy1 0; echo $xx2 $yy2 0; echo 0 0 $z2) | vaspkit | grep 123
    mv SUPERCELL.vasp POSCAR22
fi    
python3 build.py

