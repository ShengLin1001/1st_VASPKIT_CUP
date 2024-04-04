      Program read_poscar
      implicit none
      integer   i,j,k,m
      character*8  sysname_A,sysname_B
      real*8    alattice_a,AL_A(3,3),alattice_b,AL_B(3,3),alattice_C
      real*8    alatt_1,alatt_2,Rlatt_a,sta_a,Rlatt_b,sta_b
      real*8    pos_a(5000,3),pos_b(5000,3) 
      real*8    strain
      integer   ntot,natom_a,natom_b
      character(len=12) ::  form,name
      real*8   height_v,distance,delta_D,off_cen(3)
      character*255 :: line1,line2
      character*1   :: tchar
      integer   ntypeA,ntypeB,nions
      integer,allocatable :: nitypeA(:),nitypeB(:),itype(:)
      character*2,allocatable :: ele_nameA(:),ele_nameB(:)

      write(*,*)"==================================================
     &====================="
      write(*,*) "Author name: Jianwei Wang"
      write(*,*) "Email: Jianwei_wang@uestc.edu.cn"
      write(*,*) "Affiliation: University of Electronic Science and 
     &Technology of China"
      write(*,*) "Usage description:"
      write(*,*) "1.compile: ifort -o heterjunction_building.x
     &heterojunction_building.f" 
      write(*,*) "2. prepare SUPERCELL_A.vasp and SUPERCELL_B.vasp
     & as input files"
      write(*,*) "3. run heterojunction_buliding.x 
     &and the output file is POSCAR"
      write(*,*) "version: 1.0"
      write(*,*) "Last Modified Time: 2019.07.30"
      write(*,*)"==================================================
     &====================="

    
      write(*,*) 
      write(*,*) "input the height of vallcum:"
      read(*,*)  height_v
      write(*,*) "input the distance between nearby layers"
      read(*,*)  distance
      open(unit=17,file="SUPERCELL_A.vasp")
      rewind(17)
      read(17,*) sysname_A
      write(*,*) sysname_A
      read(17,*) alattice_a
c     write(*,*) alattice_a
      do i=1,3
         read(17,*) AL_A(1,i),AL_A(2,i),AL_A(3,i)
c        write(*,*)  AL_A(1,i),AL_A(2,i),AL_A(3,i)
      enddo 
      alatt_1=alattice_a*sqrt(AL_A(1,1)**2+AL_A(2,1)**2+AL_A(3,1)**2)
      write(*,*) alatt_1
      alatt_2=alattice_a*sqrt(AL_A(1,2)**2+AL_A(2,2)**2+AL_A(3,2)**2)
      write(*,*) alatt_2
      alattice_a = alatt_1
      read(17,'(A)') line1
      read(17,'(A)') line2
      call num_elements(line2,255,ntypeA) 
      write(*,'(3X,A,I5,A)')'Found ',ntypeA,' atom type(s)'
       allocate(nitypeA(ntypeA),ele_nameA(ntypeA))
             read(line1,*)(ele_nameA(i),i=1,ntypeA)
             read(line2,*)(nitypeA(i),i=1,ntypeA)
          write(*,'(3X,A)')"Number of atom for each type:"
          write(*,'(A)')"     Type   No."
          do i=1,ntypeA
             write(*,'(5X,I3,I6)')i,nitypeA(i)
          enddo
           nions=0
           do i=1,ntypeA
              nions=nions+nitypeA(i)
           enddo
          write(*,'(5X,A,I5)')"Total No. atom:",nions 
       natom_a = nions
c      read(17,*) 
c      read(17,*)  natom_a
       read(17,*)
      do i=1,natom_a
      read(17,*) pos_a(i,1),pos_a(i,2),pos_a(i,3)      
      end do
      write(*,*) "A off_centre = (x,y,z) as new origine"
      read(*,*) off_cen(1),off_cen(2),off_cen(3)      
      do i=1,natom_a
      pos_a(i,1)=pos_a(i,1)-off_cen(1)
      pos_a(i,2)=pos_a(i,2)-off_cen(2)
      pos_a(i,3)=pos_a(i,3)-off_cen(3)
      end do
      do i=1,natom_a
      pos_a(i,3)=pos_a(i,3)*AL_A(3,3)
      end do
      write(*,*) max(1,3,5,7,9)

      open(unit=18,file="SUPERCELL_B.vasp")
      rewind(18)
      read(18,*) sysname_B
      write(*,*) sysname_B
      read(18,*) alattice_b
      do i=1,3
         read(18,*) AL_B(1,i),AL_B(2,i),AL_B(3,i)
      enddo     
      alatt_1=alattice_b*sqrt(AL_B(1,1)**2+AL_B(2,1)**2+AL_B(3,1)**2)
      write(*,*) alatt_1
      alatt_2=alattice_b*sqrt(AL_B(1,2)**2+AL_B(2,2)**2+AL_B(3,2)**2)
      write(*,*) alatt_2
      alattice_b = alatt_1

      read(18,'(A)') line1
      read(18,'(A)') line2
      call num_elements(line2,255,ntypeB)
      write(*,'(3X,A,I5,A)')'Found ',ntypeB,' atom type(s)'
       allocate(nitypeB(ntypeB),ele_nameB(ntypeB))
             read(line1,*)(ele_nameB(i),i=1,ntypeB)
             read(line2,*)(nitypeB(i),i=1,ntypeB)
          write(*,'(3X,A)')"Number of atom for each type:"
          write(*,'(A)')"     Type   No."
          do i=1,ntypeB
             write(*,'(5X,I3,I6)')i,nitypeB(i)
          enddo
           nions=0
           do i=1,ntypeB
              nions=nions+nitypeB(i)
           enddo
          write(*,'(5X,A,I5)')"Total No. atom:",nions
       natom_b=nions     
             
C      read(18,*)
C      read(18,*)  natom_b1,natom_b2
      read(18,*)
      do i=1,natom_b
      read(18,*) pos_b(i,1),pos_b(i,2),pos_b(i,3)
      end do
      write(*,*) "B off_centre = (x,y,z) as new origine"
      read(*,*) off_cen(1),off_cen(2),off_cen(3)
      do i=1,natom_b
      pos_b(i,1)=pos_b(i,1)-off_cen(1)
      pos_b(i,2)=pos_b(i,2)-off_cen(2)
      pos_b(i,3)=pos_b(i,3)-off_cen(3) 
      end do    
      do i=1,natom_b
      pos_b(i,3)=pos_b(i,3)*AL_B(3,3)
      end do
      

      write(*,*)'lattice mismatch is:', alattice_a/alattice_b-1.000
      write(*,*) "combination of A and B:", sysname_A,"_",sysname_B
      alatt_1 = (alattice_a+ alattice_b)/2.0
      
      ntot=0
      alattice_C=height_v+distance               !+maxval(pos_b(1:natom_b,3))
      delta_D=maxval(pos_a(1:natom_a,3))-minval(pos_a(1:natom_a,3))
      write(*,*) 'delta A distance is :', delta_D
      alattice_C=alattice_C+delta_D
      delta_D=maxval(pos_b(1:natom_b,3))-minval(pos_b(1:natom_b,3))
      write(*,*) 'delta B distance is :', delta_D
      alattice_C=alattice_C+delta_D

      delta_D=distance-minval(pos_a(1:natom_a,3))
     &                     +maxval(pos_b(1:natom_b,3))

      open(unit=20,file="POSCAR")
      rewind(20)
      write(20,'(A6,A6)') sysname_a,sysname_b
      write(20,'(F12.8)') 1.000000
      DO I=1,2
      write(20,102) (AL_A(1,I)+AL_B(1,I))*0.50,
     &   (AL_A(2,I)+AL_B(2,I))*0.50,
     &   (AL_A(3,I)+AL_B(3,I))*0.50
      END DO

      write(20,102) 0.00,0.00,alattice_C
      write(20,'(5(A4,2X))')   (ele_nameA(i),i=1,ntypeA),
     &                         (ele_nameB(i),i=1,ntypeB)
      write(20,'(5(I4,2X))')   (nitypeA(i),i=1,ntypeA),
     &                         (nitypeB(i),i=1,ntypeB)

c      write(20,'(A6,A6,A6)') 'C','Mo','S'
c      write(20,'(I4,I4,I4)') natom_a,natom_b
      write(20,'(A8)')'Direct'
      do i=1,natom_a
c      pos_a(i,3)=3.0+ maxval(pos_b(1:natom_b,3))
      write(20,102) pos_a(i,1),pos_a(i,2),
     &                  (pos_a(i,3)+delta_D)/alattice_C
      end do
      do i=1,natom_b
      write(20,102) pos_b(i,1),pos_b(i,2),pos_b(i,3)/alattice_C
      end do


 
 101  format(F16.10,F16.10,I4,I4)
 102  format(F16.12,F16.12,F16.12)
      END 

      subroutine num_elements(charin,charl,nel)
         IMPLICIT NONE
         INTEGER :: nel,i,charl
         CHARACTER(charl):: charin
         LOGICAL ::last_int=.FALSE.
         nel =0
         do i=1,charl
            if (charin(i:i)>='0' .and. charin(i:i)<='9' .or.
     &   charin(i:i).eq.'.') then
               last_int=.TRUE.
            else if (last_int.eqv..TRUE.) then
               nel=nel+1
               last_int=.FALSE.
           endif
        end do
        end subroutine num_elements

