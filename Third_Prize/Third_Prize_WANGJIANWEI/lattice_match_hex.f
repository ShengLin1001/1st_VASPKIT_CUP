      Program read_poscar
      implicit none
      integer   i,j,k,m
      character*8  sysname_A,sysname_B
      real*8    alattice_a,AL_A(3,3),alattice_b,AL_B(3,3)
      real*8    alatt_1,alatt_2,Rlatt_a,sta_a,Rlatt_b,sta_b
      real*8    strain
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
      write(*,*) "2. prepare input files: 02_Graphene_POSCAR and
     & 02_MoS2_POSCAR"
      write(*,*) "3. run lattice_match_rect.x and the output files are
     & TRAN_A1,TRAN_B1,TRAN_A2,TRAN_B2,...."
      write(*,*) "version: 1.0"
      write(*,*) "Last Modified Time: 2019.07.30"
      write(*,*)"==================================================
     &====================="




      open(unit=17,file="02_Graphene_POSCAR")
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


      open(unit=18,file="02_MoS2_POSCAR")
      rewind(18)
      read(18,*) sysname_B
      write(*,*) sysname_B
      read(18,*) alattice_b
      do i=1,3
         read(18,*) AL_B(1,I),AL_B(2,i),AL_B(3,i)
      enddo     
      alatt_1=alattice_b*sqrt(AL_B(1,1)**2+AL_B(2,1)**2+AL_B(3,1)**2)
      write(*,*) alatt_1
      alatt_2=alattice_b*sqrt(AL_B(1,2)**2+AL_B(2,2)**2+AL_B(3,2)**2)
      write(*,*) alatt_2
      alattice_b = alatt_1
      

      write(*,*)'lattice mismatch is:', alattice_a/alattice_b-1.000

      ntot=0

      do i=1,3
      do j=0,3
      Rlatt_a=sqrt( real(i*i-i*j+j*j) )*alattice_a
      Sta_a=ACOS((real(i)-0.5*real(j))
     &         /sqrt(real(i*i-i*j+j*j)))/dacos(0.d0)*90.0
c      write(*,101) Rlatt_a,Sta_a,i,j
         
      do k=1,3
      do m=0,3
      Rlatt_b=sqrt( real(k*k-k*m+m*m) )*alattice_b
      Sta_b=ACOS((real(k)-0.5*real(m))
     &         /sqrt(real(k*k-k*m+k*m)))/dacos(0.d0)*90.0
c      write(*,101) Rlatt_b,Sta_b,k,m

      strain=abs(Rlatt_a/Rlatt_b-1.0)
c      write(*,*)  'strain =:',strain 
      if(strain <= 0.05 ) then
          write(*,*) strain, i,j,k,m
           ntot=ntot+1
           if(ntot>0  .and. ntot< 10) then
           write(form,'(i1)') ntot
           else if (ntot >= 10 .and. ntot < 100) then
           write(form,'(i2)') ntot
           else if (ntot >= 100 .and. ntot < 1000) then
           write(form,'(i3)') ntot
           end if     
           write(name,*) "TRAN_A",trim(form)
c           write(*,*) 'name = ', name
c           write(*,*) 'ntot is', ntot
         open(unit=19,file=name) 
         rewind(19)
         write(19,'(A24)') '# TRANSFORMATION  MATRIX' 
         write(19,'(I4,I4,I4)') i,j,0
         write(19,'(I4,I4,I4)') -j,i-j,0
         write(19,'(I4,I4,I4)') 0,0,1
         write(name,*) "TRAN_B",trim(form)
         open(unit=20,file=name)
         rewind(20)
         write(20,'(A24)') '# TRANSFORMATION  MATRIX'
         write(20,'(I4,I4,I4)') k,m,0
         write(20,'(I4,I4,I4)') -m,k-m,0
         write(20,'(I4,I4,I4)') 0,0,1         

        write(*,'(I4,I4,I4,A1,I4,I4,I4)') i,j,0,'|',k,m,0
        write(*,'(I4,I4,I4,A1,I4,I4,I4)') -j,i-j,0,'|',-m,k-m,0
        write(*,'(I4,I4,I4,A1,I4,I4,I4)') 0,0,1,'|',0,0,1

       end if
       end do
       end do
c      write(*,*) 'ntot is', ntot
      end do
      end do       
      
 101  format(F16.10,F16.10,I4,I4)


      END 
