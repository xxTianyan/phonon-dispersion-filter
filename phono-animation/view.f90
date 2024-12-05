!----------------------------Write by Yixuan Xue----------------------------
!----------------------------------2023.3.3---------------------------------
program main
!use read_band 
implicit none

integer :: i
integer :: j
integer :: z
integer :: n
integer :: k
integer :: q
integer :: numfile
integer :: filestati
integer :: filestatj
integer :: jdx, jdy, jdz
integer :: cot, cotsq, cotsqs, cottime, cotp      !Count
integer :: nqpoint  !Number of phonon points
integer :: npath    !Number of phonon paths
integer :: numf     !Number of frequencies
integer :: natom    !Number of atoms in crystal cell
integer :: nsatom   !Number of atoms in square cell
integer :: nstatom  !Number of atoms in square structure
integer :: nmotion  !Number of motions
integer :: numcell(3)  !Supercell_matrix in the lattice direction
integer :: numcellsquare(3)  !Supercell_matrix with a square cell
integer :: step     !Number of frames in a cycle
integer :: extendnatom !Atom number of extended cell
integer :: nmqpot
integer :: nmf
integer, allocatable :: labelexcd(:)
integer, allocatable :: labelscd(:)
integer, allocatable :: labelsqcd(:)
integer, allocatable :: nindex(:)

real(kind=8)    :: reciprocal_lattice(3,3)  !Reciprocal_lattice orientation
real(kind=8)    :: lattice(3,3)  !Lattice orientation
real(kind=8)    :: amplitude     !amplitude
real(kind=8)    :: boxcell(3)    !Cell size
real(kind=8)    :: Pi            !Pi
real(kind=8)    :: pro_reciprocal_lattice(3,3) 
real(kind=8)    :: time
real(kind=8), allocatable :: coordinates(:,:)  !Fractional coordinates of atoms
real(kind=8), allocatable :: coordinates1(:,:) 
real(kind=8), allocatable :: nq_position(:,:)  !Direction of cell movement
real(kind=8), allocatable :: nfrequency(:,:)   !Frequency of cell movement
real(kind=8), allocatable :: neigenvector(:,:,:,:,:) !Fractional coordinates of atoms
real(kind=8), allocatable :: extendcoord(:,:)  !Coordinate of extended cell
real(kind=8), allocatable :: scoord(:,:)       !Coordinate of square cell
real(kind=8), allocatable :: extendscoord(:,:) !Coordinate of extended square cell
real(kind=8), allocatable :: extendsqcoord(:,:) !Coordinate of extended square structure
real(kind=8), allocatable :: q_p(:,:) 
real(kind=8), allocatable :: f(:,:)
real(kind=8), allocatable :: period(:,:)
real(kind=8), allocatable :: dt(:,:)
real(kind=8), allocatable :: nmq_position(:,:)
real(kind=8), allocatable :: nmeigenvector(:,:,:,:,:)
real(kind=8), allocatable :: nmfrequency(:,:)
real(kind=8)    :: e_x
real(kind=8)    :: e_y
real(kind=8)    :: e_z
real(kind=8)    :: u_x
real(kind=8)    :: u_y
real(kind=8)    :: u_z
real(kind=8)    :: pos_x
real(kind=8)    :: pos_y
real(kind=8)    :: pos_z

complex(kind=8) :: e_xj
complex(kind=8) :: e_yj
complex(kind=8) :: e_zj
complex(kind=8) :: ee
complex(kind=8) :: u_px
complex(kind=8) :: u_py
complex(kind=8) :: u_pz

character :: sectype
character :: str, str1, str2
character(len=100) :: filename,filenamei,filenamej
character, allocatable :: symbol(:) !Atomic species
character, allocatable :: extendsymbol(:) !Atomic species
character, allocatable :: sqsymbol(:) !Atomic species of square cell
character, allocatable :: extendsqsymbol(:) !Atomic species of square structure

logical isExist

pi = 3.1415926

open(11,file='a.inp')
read(11,*) numcell
read(11,*) nsatom
read(11,*) boxcell
read(11,*) numcellsquare
read(11,*) amplitude
read(11,*) step
read(11,*) sectype

open(10,file = 'band.yaml')
read(10,*) str, nqpoint
read(10,*) str, npath
do i = 1, npath+2
   read(10,*)
end do
read(10,*) str, str1, reciprocal_lattice(1,:)
read(10,*) str, str1, reciprocal_lattice(2,:)
read(10,*) str, str1, reciprocal_lattice(3,:)
pro_reciprocal_lattice = reciprocal_lattice*2.0*pi
read(10,*) str, natom
read(10,*)
read(10,*) str, str1, lattice(1,1), lattice(2,1), lattice(3,1)
read(10,*) str, str1, lattice(1,2), lattice(2,2), lattice(3,2)
read(10,*) str, str1, lattice(1,3), lattice(2,3), lattice(3,3)
read(10,*)
numf = natom*3
allocate(coordinates(natom,3),coordinates1(natom,3))
allocate(symbol(natom))
allocate(nq_position(nqpoint,3))
allocate(nfrequency(nqpoint,numf))
allocate(neigenvector(nqpoint,numf,natom,3,2))
do i = 1, natom
   read(10,*) str, str1, symbol(i)
   read(10,*) str, str1, coordinates(i,:)
   read(10,*)
end do
do i = 1,6
   read(10,*)
end do

do i = 1, nqpoint
   read(10,*) str, str1, str2, nq_position(i,:)
   read(10,*)
   read(10,*)
   do j = 1, numf
      read(10,*)
	  read(10,*) str, nfrequency(i,j)
	  read(10,*)
	  do z = 1, natom
	     read(10,*)
		 do n = 1, 3
		    read(10,*) str, str1, neigenvector(i,j,z,n,:)
		 end do
	  end do
   end do
   read(10,*)
end do

extendnatom=natom*numcell(1)*numcell(2)*numcell(3)
allocate(extendcoord(extendnatom,3))
allocate(extendsymbol(extendnatom))
allocate(labelexcd(extendnatom))
jdx = INT(numcell(1)*0.5)
jdy = INT(numcell(2)*0.5)
jdz = INT(numcell(3)*0.5)

!open(12,file='text.xyz')
!write(12,*) extendnatom
!write(12,*)
cot = 0
do n = 1, natom
   do i = 1, numcell(1) 
      do j = 1, numcell(2)
         do z = 1, numcell(3)	     
	        cot = cot + 1
			extendsymbol(cot) = symbol(n)
			coordinates1(n,1)=coordinates(n,1)+(i-1)-jdx
			coordinates1(n,2)=coordinates(n,2)+(j-1)-jdy
			coordinates1(n,3)=coordinates(n,3)+(z-1)-jdz
	        extendcoord(cot,1)=Dot_product(coordinates1(n,:),lattice(1,:))
			extendcoord(cot,2)=Dot_product(coordinates1(n,:),lattice(2,:))
			extendcoord(cot,3)=Dot_product(coordinates1(n,:),lattice(3,:))
			labelexcd(cot)=n
			!write(12,*) extendsymbol(cot), extendcoord(cot,:)
	     end do
	  end do
   end do
end do

!write(12,*) nsatom
!write(12,*)
allocate(scoord(nsatom,3))
allocate(sqsymbol(nsatom))
allocate(labelscd(nsatom))
cotsq = 0
do i = 1, cot
   if ((extendcoord(i,1)>=0.0 .and. extendcoord(i,1)<boxcell(1)) .and. &
      &(extendcoord(i,2)>=0.0 .and. extendcoord(i,2)<boxcell(2)) .and. &
	  &(extendcoord(i,3)>=0.0 .and. extendcoord(i,3)<boxcell(3))) then
      cotsq=cotsq+1
	  sqsymbol(cotsq)=extendsymbol(i)
	  scoord(cotsq,:)=extendcoord(i,:)
	  labelscd(cotsq)=labelexcd(i)
	  !write(12,*) sqsymbol(cotsq), scoord(cotsq,:)
   end if
end do

nstatom = nsatom*numcellsquare(1)*numcellsquare(2)*numcellsquare(3)
!write(12,*) nstatom
!write(12,*)
allocate(extendsqcoord(nstatom,3))
allocate(extendsqsymbol(nstatom))
allocate(labelsqcd(nstatom))
cotsqs = 0
do n = 1, nsatom
   do i = 1, numcellsquare(1)
      do j = 1, numcellsquare(2)
         do z = 1, numcellsquare(3)	     
		    cotsqs=cotsqs+1
		    extendsqsymbol(cotsqs)=sqsymbol(n)
			extendsqcoord(cotsqs,1)=scoord(n,1)+(i-1)*boxcell(1)
			extendsqcoord(cotsqs,2)=scoord(n,2)+(j-1)*boxcell(2)
			extendsqcoord(cotsqs,3)=scoord(n,3)+(z-1)*boxcell(3)
			labelsqcd(cotsqs)=labelscd(n)
			!write(12,*) extendsqsymbol(cotsqs), extendsqcoord(cotsqs,:)
	     end do
	  end do
   end do
end do

if (sectype=='W') then
    nmqpot = nqpoint
	nmf = numf
	allocate(nmq_position(nmqpot,3))
    allocate(nmfrequency(nmqpot,nmf))
    allocate(nmeigenvector(nmqpot,nmf,natom,3,2))
	nmfrequency = nfrequency
	nmq_position = nq_position
	nmeigenvector = neigenvector
else if (sectype=='S') then
    read(11,*) nmqpot
	allocate(nindex(nmqpot))
	read(11,*) nindex
    nmf = numf
	allocate(nmq_position(nmqpot,3))
    allocate(nmfrequency(nmqpot,nmf))
    allocate(nmeigenvector(nmqpot,nmf,natom,3,2))
	do i = 1, nmqpot
       nmq_position(i,:) = nq_position(nindex(i),:)
	   nmfrequency(i,:) = nfrequency(nindex(i),:)
	   nmeigenvector(i,:,:,:,:) = neigenvector(nindex(i),:,:,:,:)
	end do
else if (sectype=='Z') then
   read(11,*) nmqpot
   nmf = 1
   allocate(nmq_position(nmqpot,3))
   allocate(nmfrequency(nmqpot,nmf))
   allocate(nmeigenvector(nmqpot,nmf,natom,3,2))
   do i = 1, nmqpot
      read(11,*) nmq_position(i,:)
	  read(11,*) nmfrequency(i,1)
	  do j = 1, natom
	     read(11,*) nmeigenvector(i,nmf,j,:,:)
	  end do
   end do
end if

allocate(q_p(nmqpot,3))
allocate(f(nmqpot,nmf))
allocate(period(nmqpot,nmf))
allocate(dt(nmqpot,nmf))
e_x = 0.0
e_y = 0.0
e_z = 0.0
e_xj = Cmplx(0.0, 0.0)
e_yj = Cmplx(0.0, 0.0)
e_zj = Cmplx(0.0, 0.0)
ee = Cmplx(0.0, 0.0)
u_px = Cmplx(0.0, 0.0)
u_py = Cmplx(0.0, 0.0)
u_pz = Cmplx(0.0, 0.0)
u_x = 0.0
u_y = 0.0
u_z = 0.0
pos_x = 0.0
pos_y = 0.0
pos_z = 0.0
f=nmfrequency*2.0*pi
period=2.0*pi/f
dt=period/Real(step)
numfile=13

inquire(file='Results',exist=isExist)
if (isExist) then
   call system ("rd/s/q "//trim("Results"))   !For Windows
!   call system ("rm "//trim("Results"))   !For Linux
else
   call system("mkdir "//trim("Results"))
end if

do i = 1, nmqpot
   filestati=i
   write(filenamei,'(i5)') filestati
   filenamei=adjustl(filenamei)
   q_p(i,1)=Dot_product(pro_reciprocal_lattice(1,:), nmq_position(i,:))
   q_p(i,2)=Dot_product(pro_reciprocal_lattice(2,:), nmq_position(i,:))
   q_p(i,3)=Dot_product(pro_reciprocal_lattice(3,:), nmq_position(i,:))
   call system("mkdir " //trim("Results\")//trim(filenamei))
   do j = 1, nmf
      !numfile=numfile+1	  
	  filestatj=j
	  write(filenamej,'(i5)') filestatj
	  filenamej=adjustl(filenamej)
	  filename=trim("Results\")//trim(filenamei)//"\"//trim(filenamej)//'.xyz'
	  open(numfile,file=filename)
      do z = 1, step
         cottime=-1
         time=Real(z-1)*dt(i,j)
		 write(numfile,'(g0)') nstatom
         write(numfile,*) "Timestep:", z		 
	     do k = 1, nstatom
		    n=labelsqcd(k) 			
		    e_x=nmeigenvector(i,j,n,1,1)
            e_xj=Cmplx(0.0, nmeigenvector(i,j,n,1,2))
            e_y=nmeigenvector(i,j,n,2,1)
            e_yj=Cmplx(0.0, nmeigenvector(i,j,n,2,2))
            e_z=nmeigenvector(i,j,n,3,1)
            e_zj=Cmplx(0.0, nmeigenvector(i,j,n,3,2))
            ee=Exp(Cmplx(0.0,time*f(i,j)-Dot_product(extendsqcoord(k,:),q_p(i,:))))
			u_px=(e_x+e_xj)*ee
            u_py=(e_y+e_yj)*ee
            u_pz=(e_z+e_zj)*ee
            
            u_x=Real(u_px)*amplitude
            u_y=Real(u_py)*amplitude
            u_z=Real(u_pz)*amplitude
                            
            pos_x=extendsqcoord(k,1)+u_x
            pos_y=extendsqcoord(k,2)+u_y
            pos_z=extendsqcoord(k,3)+u_z
			write(numfile,*) extendsqsymbol(k), pos_x, pos_y, pos_z
         end do
      end do
   close(numfile)
   end do
end do

write(*,*) 'Good job! See you later~^-^'
end