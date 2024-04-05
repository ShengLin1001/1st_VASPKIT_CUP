#!/opt/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:16:57 2019

@author: XueFei Liu 316187631@qq.com
"""
import re
import sys
import os
from numpy import *
import numpy as np
import time


def lattice_write(file):
    lattice=[]
    poscar_tem=[]
    scale=[]
    if os.path.exists(file):
       poscar=open(file,'r')
    
    else:
        raise IOError('POSCAR does not exist！')
    for line in poscar.readlines():
                poscar_tem.append(line)
    poscar.close()
    scale=poscar_tem[1].split()
    lat=poscar_tem[2:5]
    for i in lat:
        lattice.append(i.split())
    scale_lat=[]
    scale_lat.append(scale)
    scale_lat.append(lattice)
    return scale_lat
def run_vaspkit(mna=[1,0],mnb=[0,1],poscar='POSCAR'):
        os.environ['poscar']=str(poscar)
    
        ma_x0=mna[0]
        ma_y0=mna[1]
        os.environ['ma_x0']=str(ma_x0)
        os.environ['ma_y0']=str(ma_y0)
        ma_x1=mnb[0]
        ma_y1=mnb[1]
        os.environ['ma_x1']=str(ma_x1)
        os.environ['ma_y1']=str(ma_y1)
        os.system('(echo $ma_x0 $ma_y0 0;\
                  echo $ma_x1 $ma_y1 0;echo 0 0 1)|vaspkit -task 400\
                  -file $poscar > /dev/null')
def     searchRange(nums, target):
       
        flag = 0
        lis = []
        for i in range(nums.count(target)):
            sec = flag
            flag = nums[flag:].index(target)
            lis.append(flag + sec)
            flag = lis[-1:][0] + 1
        return lis
if __name__=='__main__':
    start=time.time()
    mismatch=0.05                 
    print("*************************************************************")
    print("Author:XueFei Liu")
    print("Email:316198731@qq.com")
    print("Function:A code to build 2D hetero-junction with two POSCAR")
    print("Version:01")
    print("First release time:2019.7.11")
    print("*************************************************************")
   
    poscar_a= sys.argv[1]
    poscar_b= sys.argv[2]
    scal_lata=lattice_write(poscar_a)
    scal_latb=lattice_write(poscar_b)
    a0=[]
    b0=[]
    a1=[]
    b1=[]
    
    for i in scal_lata[1][0]:
        a0.append(float(i))
    
    for i in scal_lata[1][1]:
        b0.append(float(i))
    
    for i in scal_latb[1][0]:
        a1.append(float(i))
   
    for i in scal_latb[1][1]:
        b1.append(float(i))
             
    a2=np.sqrt(float(scal_lata[1][0][0])**2+\
               float(scal_lata[1][0][1])**2+\
               float(scal_lata[1][0][2])**2)
    b2=np.sqrt(float(scal_latb[1][0][0])**2+\
               float(scal_latb[1][0][1])**2+\
               float(scal_latb[1][0][2])**2)
    print("*************************************************************")
    print("   Input the scale times relative to old lattice,")
    print("(  Noting:You'd better increase the value from a smaller one e.g:2 3 ...")
    print("   A smaller value will run speeder and produce an ")
    print("   as smaller supercell as possible)" )
    print("*************************************************************")
    evalue=int(input("Input scale times:"))
    try:
      mismatch=float(input("Input the mismatch of two lattice:"))
    except:
        print("The input is invalid,the default of mismatch is 0.05")
        mismatch=0.05
    delta=0.3
    a0_new_f=[]
    b0_new_f=[]
    a1_new_f=[]
    b1_new_f=[]
    ab_angle_new_f=[]
    ab0_new=[]
    a0_new_v=[]
    b0_new_v=[]
    c_ab0_new=[]
    angle_ab0_new=[]
    ab1_new=[]
    a1_new_v=[]
    b1_new_v=[]
    c_ab1_new=[]
    angle_ab0_new=[]
    angle_ab0_new0=[]
    angle_ab1_new1=[]
    a0_new_f0=[]
    b0_new_f0=[]
    a1_new_f1=[]
    b1_new_f1=[]
    a0_new_f01=[]
    b0_new_f01=[]
    a1_new_f11=[]
    b1_new_f11=[]
    ab_angle_new_f=[]
    ab_angle_new_f1=[]
    ab_angle_new_f2=[]
    ab_angle_new_f_min=[]
    m_n=[]
    m_n_a0_new=[]
    m_n_b0_new=[]
    m_n_a1_new=[]
    m_n_b1_new=[]
    m_n_a0_newf=[]
    m_n_b0_newf=[]
    m_n_a1_newf=[]
    m_n_b1_newf=[]
    a0_new_f_index=[]
    a1_new_f_index=[]
    b0_new_f_index=[]
    b1_new_f_index=[]
    a0_new_v_index=[]
    a1_new_v_index=[]
    b0_new_v_index=[]
    b1_new_v_index=[]
    a0_v=round(np.linalg.norm(a0),1)
    b0_v=round(np.linalg.norm(b0),1)
    a1_v=round(np.linalg.norm(a1),1)
    b1_v=round(np.linalg.norm(b1),1)
    for m in range(-evalue,evalue):
        for n in range(-evalue,evalue):
          if m ==0 and n == 0:
              continue
          else:
            
            a0_new_f.append([round(m*a0[i] + n*b0[i],1) for i in range(3)])
            b0_new_f.append([round(m*a0[i] + n*b0[i],1) for i in range(3)])
            
           
            a1_new_f.append([round(m*a1[i] + n*b1[i],1) for i in range(3)])
            b1_new_f.append([round(m*a1[i] + n*b1[i],1) for i in range(3)])
            m_n.append([m,n])
           
    for i in range(len(a0_new_f)):
           
        for j in range(len(b0_new_f)):
             
             ab0_new.append(round(np.dot(a0_new_f[i],b0_new_f[j]),1))
             a0_new_v.append(round(np.linalg.norm(a0_new_f[i]),1))
             b0_new_v.append(round(np.linalg.norm(b0_new_f[j]),1))
             a0_new_f_index.append(a0_new_f[i])
             b0_new_f_index.append(b0_new_f[j])
                
    for i in range(len(a1_new_f)):
           
        for j in range(len(b1_new_f)):
             
             ab1_new.append(round(np.dot(a1_new_f[i],b1_new_f[j]),1))
             a1_new_v.append(round(np.linalg.norm(a1_new_f[i]),1))
             b1_new_v.append(round(np.linalg.norm(b1_new_f[j]),1)) 
             a1_new_f_index.append(a1_new_f[i])
             b1_new_f_index.append(b1_new_f[j])
    
    for i in range(len(ab0_new)):
        angle_ab0_new0.append(round(np.arccos(ab0_new[i]/(a0_new_v[i]*b0_new_v[i]))*180/np.pi,0))
    for i in range(len(ab1_new)):
        angle_ab1_new1.append(round(np.arccos(ab1_new[i]/(a1_new_v[i]*b1_new_v[i]))*180/np.pi,0))
    for i in range(len(a0_new_v)):
        for j in range(len(a1_new_v)):
            if abs((a0_new_v[i]-a1_new_v[j])/a0_new_v[i]) < mismatch and \
            abs(angle_ab0_new0[i] - angle_ab1_new1[j]) < delta and \
            abs((b0_new_v[i]-b1_new_v[j])/b0_new_v[i]) < mismatch and \
            angle_ab0_new0[i] > 30 and angle_ab0_new0[i]<150 and \
            angle_ab1_new1[j] > 30 and angle_ab1_new1[j] < 150 :
               a0_new_f0.append(round(a0_new_v[i],1))
               
               a1_new_f1.append(round(a1_new_v[j],1))
               b0_new_f0.append(round(b0_new_v[i],1))
               
               b1_new_f1.append(round(b1_new_v[j],1))
               ab_angle_new_f1.append(round(angle_ab0_new0[i],0))
               ab_angle_new_f2.append(round(angle_ab1_new1[j],0))
               a0_new_v_index.append(a0_new_f_index[i])
              
               a1_new_v_index.append(a1_new_f_index[j])
               b0_new_v_index.append(b0_new_f_index[i])
               
               b1_new_v_index.append(b1_new_f_index[j])
    a0_new_f0_dup=[]
    b0_new_f0_dup=[]
    a1_new_f1_dup=[]
    b1_new_f1_dup=[]
    if len(a0_new_f0) !=0:            
       min_lat=min(a0_new_f0)
    else:
       print("Please increase scale times and rerun this code!")
       sys.exit(0)
    a0_new_f0_min=a0_new_f0[:]
    for i in range(len(a0_new_f0_min)):
        if abs(a0_new_f0_min[i]-min_lat)<0.3:
           a0_new_f0_min[i]=min_lat 
    min_index=searchRange(a0_new_f0_min,min_lat)
    for i in min_index:
        m_n_a0_new.append(m_n[a0_new_f.index(a0_new_v_index[i])])
        m_n_b0_new.append(m_n[a0_new_f.index(b0_new_v_index[i])])
        
        m_n_a1_new.append(m_n[a1_new_f.index(a1_new_v_index[i])])
        m_n_b1_new.append(m_n[b1_new_f.index(b1_new_v_index[i])])
        ab_angle_new_f_min.append(ab_angle_new_f1[i])
        a0_new_f0_dup.append(a0_new_f0[i])
        b0_new_f0_dup.append(b0_new_f0[i])
        a1_new_f1_dup.append(a1_new_f1[i])
        b1_new_f1_dup.append(b1_new_f1[i])
    
    print("*****************************************************")
    print("The possible m-n couples ：")
    print("m_n_a1:{}".format(m_n_a0_new))
    print("m_n_b1:{}".format(m_n_b0_new))
    print("m_n_a2:{}".format(m_n_a1_new))
    print("m_n_b2:{}".format(m_n_b1_new))
    print("*****************************************************")
    print("The new lattice length: ")
    print("New_a1:{}".format(a0_new_f0_dup))
    print("New_b1:{}".format(b0_new_f0_dup))
    print("New_a2:{}".format(a1_new_f1_dup))
    print("New_b2:{}".format(b1_new_f1_dup))
    print("*****************************************************")
    print("The possible lat-angle of a-b in current scale times:")    
    print("*****************************************************")
   
    print(ab_angle_new_f_min)
    index_angle=[]
    if len(ab_angle_new_f_min) !=0:
       
       angles=input("Input angles you want to produce new POSCAR:")
       
    else:
       print("*********************************************************************")
       print("No m_n couple found! Please increase scale times and rerun this code!")
       print("*********************************************************************")
       sys.exit(0)
    for angle in angles.split():
      index_angle.append(searchRange(ab_angle_new_f_min,round(float(angle),0)))
    for index in range(len(index_angle)):
        
      os.environ['index']=str(index)
      for i in index_angle[index]:
        os.environ['i']=str(i)
        
        run_vaspkit(mna=m_n_a0_new[i],mnb=m_n_b0_new[i],poscar=poscar_a)
        if os.path.exists("SUPERCELL.vasp"):
           os.system('mv SUPERCELL.vasp poscar_a_$index-$i')
        
        run_vaspkit(mna=m_n_a1_new[i],mnb=m_n_b1_new[i],poscar=poscar_b)
        if os.path.exists("SUPERCELL.vasp"):
           os.system('mv SUPERCELL.vasp poscar_b_$index-$i')
           
       
      print("Say thanks to vaspkit when you use this code !")   
      for i in index_angle[index]:
        
        if os.path.exists("poscar_a_"+str(index)+"-"+str(i)) and\
        os.path.exists("poscar_b_"+str(index)+"-"+str(i)):
           os.environ['i']=str(i)
           
           os.system('python heterojunction.py poscar_a_$index-$i poscar_b_$index-$i')
           os.system('mv POSCAR_NEW hetero-POSCAR_NEW_$index-$i')
           os.system('echo The new hetero-POSCAR is written in \
                  hetero-POSCAR_NEW_$index-$i')
           os.system('rm poscar_*')
           break
           
      print("The smallest hetero-POSCAR has been prouduced!")
      hetero=0
      for i in index_angle[index]:     
        if os.path.exists("hetero-POSCAR_NEW_"+str(index)+"-"+str(i)):
           hetero=1 
           break
    
      if hetero == 0:        
        print("No heterojunction is produced! please input a larger scale time !")
        
    end=time.time()
    print("The total time of this calculation is: {:d} seconds".format(int(end-start)))
