      Program read_poscar
      implicit none
      integer   i,j,k,m
      character*16  sysname_A,sysname_B
      real*8    alattice_a,AL_A(3,3),alattice_b,AL_B(3,3)
      real*8    alatt_A_1,alatt_A_2,alatt_B_1,alatt_B_2
      real*8    Rlatt_a1,Rlatt_a2,Rlatt_b1,Rlatt_b2
      real*8    strain11,strain12,strain21,strain22
      integer   ntot
      character(len=12) ::  form,name


      write(*,*)"==================================================
     &====================="
      write(*,*) "Author name: Jianwei Wang"
      write(*,*) "Email: Jianwei_wang@uestc.edu.cn"
      write(*,*) "Affiliation: University of Electronic Science and
     &Technology of China"
      write(*,*) "Usage description:"
      write(*,*) "1. compile: ifort -o lattice_match_rect.x
     &lattice_match_rect.f"
      write(*,*) "2. prepare input files: 03_black_P_POSCAR and
     & 03_Graphene_POSCAR_rect"
      write(*,*) "3. run lattice_match_rect.x and the output files are
     & TRAN_A1,TRAN_B1,TRAN_A2,TRAN_B2,...."
      write(*,*) "version: 1.0"
      write(*,*) "Last Modified Time: 2019.07.30"
      write(*,*)"==================================================
     &====================="


      open(unit=17,file="03_black_P_POSCAR")
      rewind(17)
      read(17,*) sysname_A
      write(*,*) sysname_A
      read(17,*) alattice_a
      write(*,*) alattice_a
      do i=1,3
         read(17,*) AL_A(1,i),AL_A(2,i),AL_A(3,i)
        write(*,*)  AL_A(1,i),AL_A(2,i),AL_A(3,i)
      enddo 
      alatt_A_1=alattice_a*sqrt(AL_A(1,1)**2+AL_A(2,1)**2+AL_A(3,1)**2)
      write(*,*) alatt_A_1
      alatt_A_2=alattice_a*sqrt(AL_A(1,2)**2+AL_A(2,2)**2+AL_A(3,2)**2)
      write(*,*) alatt_A_2
c      alattice_a = alatt_A_1


      open(unit=18,file="03_Graphene_POSCAR_rect")
      rewind(18)
      read(18,*) sysname_B
      write(*,*) sysname_B
      read(18,*) alattice_b
      do i=1,3
         read(18,*) AL_B(1,I),AL_B(2,i),AL_B(3,i)
      enddo     
      alatt_B_1=alattice_b*sqrt(AL_B(1,1)**2+AL_B(2,1)**2+AL_B(3,1)**2)
      write(*,*) alatt_B_1
      alatt_B_2=alattice_b*sqrt(AL_B(1,2)**2+AL_B(2,2)**2+AL_B(3,2)**2)
      write(*,*) alatt_B_2
c      alattice_b = alatt_B_1
      

c      write(*,*)'lattice mismatch is:', alattice_a/alattice_b-1.000

      ntot=0

      do i=1,10
      do j=1,10
      Rlatt_a1=real(i)*alatt_A_1
      Rlatt_a2=real(j)*alatt_A_2
c      write(*,101) Rlatt_a,Sta_a,i,j
         
      do k=1,10
      do m=1,10
      Rlatt_b1=real(k)*alatt_B_1
      Rlatt_b2=real(m)*alatt_B_2

c      write(*,101) Rlatt_b,Sta_b,k,m

      strain11=abs(Rlatt_a1/Rlatt_b1-1.0)
      strain22=abs(Rlatt_a2/Rlatt_b2-1.0)
      strain12=abs(Rlatt_a1/Rlatt_b2-1.0)
      strain21=abs(Rlatt_a2/Rlatt_b1-1.0)

c      write(*,*)  'strain =:',strain 
      if( (strain11 <=0.05) .and. (strain22 <=0.05) ) then
          write(*,*) strain11,strain22, i,j,k,m
           ntot=ntot+1
           if(ntot>0  .and. ntot< 10) then
           write(form,'(i1)') ntot
           else if (ntot >= 10 .and. ntot < 100) then
           write(form,'(i2)') ntot
           else if (ntot >= 100 .and. ntot < 1000) then
           write(form,'(i3)') ntot
           end if     
           write(name,*) "TRAN_A",trim(form)
           write(*,*) 'name = ', name
           write(*,*) 'ntot is', ntot
         open(unit=19,file=name) 
         rewind(19)
         write(19,'(A24)') '# TRANSFORMATION  MATRIX' 
         write(19,'(I4,I4,I4)') i,0,0
         write(19,'(I4,I4,I4)') 0,j,0
         write(19,'(I4,I4,I4)') 0,0,1
         write(name,*) "TRAN_B",trim(form)
         open(unit=20,file=name)
         rewind(20)
         write(20,'(A24)') '# TRANSFORMATION  MATRIX'
         write(20,'(I4,I4,I4)') k,0,0
         write(20,'(I4,I4,I4)') 0,m,0
         write(20,'(I4,I4,I4)') 0,0,1         

        write(*,'(I4,I4,I4,A1,I4,I4,I4)') i,0,0,'|',k,0,0
        write(*,'(I4,I4,I4,A1,I4,I4,I4)') 0,j,0,'|',0,m,0
        write(*,'(I4,I4,I4,A1,I4,I4,I4)') 0,0,1,'|',0,0,1

       end if
c       end do
c       end do
c      write(*,*) 'ntot is', ntot
c      end do
c      end do       
      
       if( (strain12 <=0.05) .and. (strain21 <=0.05) ) then
          write(*,*) strain11,strain22, i,j,k,m
           ntot=ntot+1
           if(ntot>0  .and. ntot< 10) then
           write(form,'(i1)') ntot
           else if (ntot >= 10 .and. ntot < 100) then
           write(form,'(i2)') ntot
           else if (ntot >= 100 .and. ntot < 1000) then
           write(form,'(i3)') ntot
           end if
           write(name,*) "TRAN_A",trim(form)
           write(*,*) 'name = ', name
           write(*,*) 'ntot is', ntot
         open(unit=19,file=name)
         rewind(19)
         write(19,'(A24)') '# TRANSFORMATION  MATRIX'
         write(19,'(I4,I4,I4)') i,0,0
         write(19,'(I4,I4,I4)') 0,j,0
         write(19,'(I4,I4,I4)') 0,0,1
         write(name,*) "TRAN_B",trim(form)
         open(unit=20,file=name)
         rewind(20)
         write(20,'(A24)') '# TRANSFORMATION  MATRIX'
         write(20,'(I4,I4,I4)') k,0,0
         write(20,'(I4,I4,I4)') 0,m,0
         write(20,'(I4,I4,I4)') 0,0,1

          write(*,'(I4,I4,I4,A1,I4,I4,I4)') i,0,0,'|',k,0,0
        write(*,'(I4,I4,I4,A1,I4,I4,I4)') 0,j,0,'|',0,m,0
        write(*,'(I4,I4,I4,A1,I4,I4,I4)') 0,0,1,'|',0,0,1

       end if
       end do
       end do
c      write(*,*) 'ntot is', ntot
       end do
       end do


     
 101  format(F16.10,F16.10,I4,I4)


      END 
