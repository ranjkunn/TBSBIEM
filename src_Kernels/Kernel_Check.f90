! MIT License

! Copyright (c) 30-04-2025 Y-Sharath-Chandra-Mouli, Ranjith-Kunath, and Avinash-Gupta

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!###################################################################################################################################
!
! Global Parameters:
! SP  = Single precision
! DP  = Double precision
! PII = The constant pi
! E6  = 10E6

! This saction documents the main code representing each variable and its use:
!  wave1 -  set up wave nos. (k -x1 direction)
!  wave3 -  set up wave nos. (k -x3 direction)
!  tu    -  T in put for calculation of upper halfsapce inplane kernels
!  tuu   -  T in put for calculation of upper halfsapce antiplane kernels
!  td    -  T in put for calculation of lower halfsapce inplane kernels
!  tdd   -  T in put for calculation of lower halfsapce antiplane kernels
!  m11u,m11d,m12u,m12d,m22u,m22d -  kernel data variables
!  qc    -  Normalized Wave number i.e. dsqrt(wave1(i)**2+wave3(j)**2)
!  slip1 -  Slip in x1 direction
!  slip2 -  Separation in x2 direction
!  slip3 -  Slip in x3 direction
!  slip  -  Slip magnitude i.e. dsqrt(slip1(i,j)**2 + slip3(i,j)**2)
!  c11u,c11d,c121u,c122u,c121d,c122d,c22u,c22d,c33u,c33d - Terms in convolution matrix
!  mus   -  Distrubution of static fraction coefficient
!  mur   -  Distrubution of dynamic fraction coefficient
!  Tau1Old,Tau2Old,Tau3Old  -  Tau's in previous iteration of Tn Convergence
!  velo1, velo2,velo3 - slip velocity distribution at given time increment
!  TT0, TS0 - Initial Distribution of Tau1 and Tau3 at 
!  Tn    - Normal stress distribution in a given time increment
!  Tn_err - Distrobution of Tn error in a given Tn convergence itteration
!  Tf    - Distribution of shear strengths at a given time increment
!  taumag-  Distribution of shear stresses at a given time increment
!  const1,const2, const3 - Temparary constants for slip rate calculation
!  conv11,conv12,conv13,conv21,conv22,conv23,conv31,conv32,conv33, conv1,conv2,conv3 - Convolution terms
!  tauf1, tauf2,tauf3 - Tau- Tnu_0 distribution at a given time increment
!  tau1, tau3,tau2 - Spatio-temparal distributions of Tau1, Tau2, and Tau3
!  TAU1F_D , TAU2F_D, TAU3F_D, CONV1_D, CONV2_D, CONV3_D - Variables for CuFFT device to host transfer
!  tauf1trans, tauf2trans, tauf3trans, conv1trans, conv2trans, conv3trans -  transpose of obtained CuFFT and inverse values
!  PLAN  - Variable to create a FFT plan for CuFFT
!  gamma - A non-dimensional time step to discretize the kernels
!  nu    - Poisson's ratio
!  time  - Time
!  betaa - Courant parameter
!  x1, x3, i, j - Spacial variables for loops
!  c111u,c111d,c1121u,c1121d,c222u,c222d,c333u,c333d - Temporary Convolution values
!  nele1 - Number of discretization points in x1 - direction 
!  nele3 - Number of discretization points in x3 - direction 
!  itter - itteration number in Tn convergence loop
!  Nitter- Maximum number of iutteration for Tn convergence loop
!  TPV_No- TPV Problem number
!  ntime - Total number time increments
!  itime - current time increment
!  ktime - temp time variable for convolution histories
!  n_pts1,n_pts3 - Number of Kernels points in X1 and X3 directions respectively
!  START_TIME, END_TIME, ELAPSED_TIME - Start time, end time and time ellpsed druing the computation
!  Tend  - Total time of computtaion
!  factor1, factor2 - Some temperary variable for convolutions
!  L1    - The total length of planer interface in X1 direction 
!  L3    - The total length of planer interface in X3 direction 
!  L1rpt - The length of slip-weakening interface size in X1 direction
!  L3rpt - The length of slip-weakening interface size in X3 direction
!  dt    - Change in time increment 
!  T0bg  - Shear stress outside of the nucleation patch (background shear stress) in (MPa).
!  Tn0   - Initial Normal stress in (MPa)
!  deltac- Critical slip-weakening distance 
!  dx    - Spatial discretization parameter
!  T0nu  - Nucleation shear stress at the center in (MPa) for TPV3, TPV5, TPV6, and TPV7.
!  T0nuL - Nucleation shear stress at the left in (MPa) for TPV5. Use T0bg for TPV3, TPV6 and TPV7.
!  T0nuR - Nucleation shear stress at the right in (MPa) for TPV5. Use T0bg for TPV3, TPV6 and TPV7.
!  Cs    - The shear wave speed of lower half-space
!  csratio-csratio of lower and top half-space
!  Csm   - Dilation wave speed of lower half-space
!  Cd    - Dilation wave speed of top half-space
!  Mu    - shear modulus of top half-space
!  Mum   - shear modulus of lower half-space
!  Muratio-shear modulus ratio of lower and top half-space
!  cdratio-cdratio of lower and top half-space
!  Cdm   - Dilation wave speed of top half-space
!  Ro    - The density of top half-space
!  Roratio-density ratio of lower and top half-space
!  Rom   - The density of lower half-space
!  css   - Cs*Cs
!  cssm  - Csm*Csm
!  eps   - Tolerance for Tn convergence norm
!  mus0  - static fraction coefficient
!  mur0  - dynamic fraction coefficient
!  outinterval1 -  Time intervel-1 for output of resutls
!  outinterval2 - Time intervel-2 for output of resutls
!  outpos3  - Spatial position-3 for output of resutls
!  outpos2  - Spatial position-2 for output of resutls
!  outpos1  - Spatial position-1 for output of resutls
!  Pre_Kernels -  Flag for precalculated kernels '0' if needed to calculate, '1' if Kernels pre calculated
!  filename - Temp variable for file names
!  word160 - Temp variable to read input file
!  pb_name - TPV problem selection variable
!  Kernel_Dir - Directes where the Kernels are stored.
!###################################################################################################################################
! Read me before running the code
! Create folders ./Kernels and ./dats
! Before running the code
! main program code begins
PROGRAM Kernel_Check
   USE Kernel

   IMPLICIT NONE   
   INTEGER, PARAMETER :: SP = KIND(1.0)
   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   real(DP),parameter :: PII = 4.d0 * datan(1.0D0)
   real(DP),parameter :: E6 = 1000000.D0
        
   ! Step 1.1: Decleration of ALLOCATABLE VARIABLES
   REAL(DP),  ALLOCATABLE, DIMENSION(:)  :: tu,td,tuu,tdd, m11u,m11d,m12u,m12d,m22u,m22d !(wave1 -x1 direction)  (wave3 -x3 direction)
   ! REAL(DP),  ALLOCATABLE, DIMENSION(:,:):: qc,slip,slip1,slip2,slip3,c11u,c11d,c121u,c122u,c121d,c122d,c22u,c22d,c33u,c33d
   REAL(DP),  ALLOCATABLE, DIMENSION(:,:):: mus,mur, station_pt
   REAL(DP),  ALLOCATABLE, DIMENSION(:,:)::T0nu
   ! REAL(DP),  ALLOCATABLE, DIMENSION(:,:)::T0nu, Tau1Old,Tau2Old,Tau3Old 
   ! Real(DP),  ALLOCATABLE, DIMENSION(:,:):: velo1, velo2,velo3, TT0, TS0, Tn, Tn_err, Tf, taumag, const1,const2, const3
   ! Complex(DP),ALLOCATABLE, DIMENSION(:,:):: conv11,conv12,conv13,conv21,conv22,conv23,conv31,conv32,conv33, conv1,conv2,conv3
   ! Complex(DP),ALLOCATABLE, DIMENSION(:,:):: tauf1, tauf2,tauf3 ! taulf(x1,x3)
   ! Complex(DP),ALLOCATABLE, DIMENSION(:,:,:):: tau1, tau3,tau2   !  tau(x1,x3,t)
   
   ! STEP 1.2: Decleration of  ALLOCATABLES OF VARIABLES FOR cuFFT
   ! Complex(DP), DEVICE, ALLOCATABLE :: TAU1F_D(:,:), CONV1_D(:,:)
   ! Complex(DP), DEVICE, ALLOCATABLE :: TAU2F_D(:,:), CONV2_D(:,:)
   ! Complex(DP), DEVICE, ALLOCATABLE :: TAU3F_D(:,:), CONV3_D(:,:)
   ! complex(DP), ALLOCATABLE , DIMENSION(:,:)  :: tauf1trans, tauf2trans, tauf3trans, conv1trans, conv2trans, conv3trans
   Real (DP), ALLOCATABLE, DIMENSION(:):: Asp_x, Asp_y, Asp_rad, Asp_mus0, Asp_mur0, T0nux11, T0nux31, T0nux12, T0nux32
   ! INTEGER :: PLAN
   ! Step 1.3: Decleration of variables
   REAL(DP) :: gamma, nu, time, betaa, x1, x3,Xr
   ! REAL(DP) :: c111u,c111d,c1121u,c1121d,c222u,c222d,c333u,c333d
   INTEGER :: itter,Nitter
   INTEGER :: ntime, nele1, nele3, itime, ktime, i, j, n_pts1,n_pts3
   REAL(DP) :: START_TIME, END_TIME, ELAPSED_TIME, Tend, factor1, factor2
   REAL(DP) :: L1,L3,L1rpt,L3rpt,dt,T0bg(3),Tn0,deltac,dx,T0nuL,T0nuR
   REAL(DP) :: Cs,csratio,Csm,Cd,Mu,Mum,Muratio,cdratio,Cdm,Ro,Roratio,Rom, css, cssm
   REAL(DP) :: eps, mus0,mur0
   INTEGER ::  outinterval1, outinterval2
   INTEGER :: outpos3,outpos2, outpos1
   INTEGER :: Pre_Kernels, io_status,No_Asp,No_neucles,NStations
   character (len=90) :: filename
   CHARACTER (len=160) :: word160   
   CHARACTER (len=20) :: pb_name, TPV_No,Kernel_Dir
   CHARACTER (len=1) :: Z_SYMM
   CALL CPU_TIME(START_TIME)   
   ! Step 1.5: Decleration of parameter

   OPEN(9001,FILE='./Pre_Kernel_status',STATUS='UNKNOWN')
   OPEN(9002,FILE='./Pre_Kernel_error',STATUS='UNKNOWN') 
! Read the test problem to be solved and open the input parameter file of corresponding problem
#ifdef TEST
      Write(9001,*)'Running in TEST mode'
      pb_name = 'TEST'
      flush(9001)
      OPEN(7008,FILE='./Input_files/TEST.in',STATUS='old',iostat=io_status)
#else
   OPEN(7001,FILE='./TPV_Problem.in',STATUS='OLD',iostat=io_status)
   if (io_status == 0) then
      READ(7001,*) pb_name 
      Write(9001,*) 'Solving ', pb_name
      flush(9001)
      IF (pb_name == 'TPV3') then
         OPEN(7008,FILE='./Input_files/TPV3.in',STATUS='unknown')
      ELSEIF (pb_name == 'TPV5') then
         OPEN(7008,FILE='./Input_files/TPV5.in',STATUS='unknown')         
      ELSEIF (pb_name == 'TPV6') then
         OPEN(7008,FILE='./Input_files/TPV6.in',STATUS='unknown')
      ELSEIF (pb_name == 'TPV7') then
         OPEN(7008,FILE='./Input_files/TPV7.in',STATUS='unknown')     
     else
        write(9001,*)  'Input Error'
        write(9002,*)  'Specify the problem to sovle: TPVx  ! where x can be 3,5,6, or 7'
        CLOSE(9001)
        CLOSE(9002)
        stop
     end if
   end if
   CLOSE(7001)
#endif
      
   ! ! Step 1.4: Read the problem input perameters file
      READ(7008,*) word160 ! Line 1
      READ(7008,*) word160 ! Line 2
      READ(7008,*) TPV_No  ! Line 3
      READ(7008,*) word160 ! Line 4
      READ(7008,*) word160 ! Line 5 
      READ(7008,*) L1rpt   ! Line 6
      READ(7008,*) word160 ! Line 7
      READ(7008,*) L3rpt   ! Line 8
      READ(7008,*) word160 ! Line 9  
      READ(7008,*) Xr      ! Line 10     
      READ(7008,*) word160 ! Line 11
      READ(7008,*) Z_SYMM    ! Line 12
      L1 = Xr * L1rpt
      if (Z_SYMM .eq. 'T') THEN     ! Check for Z_Symmetry
         L3rpt = 2.d0 * L3rpt
         L3 = Xr * L3rpt   
      else 
         L3 = Xr * L3rpt   
      endif
      READ(7008,*) word160 ! Line 13
      READ(7008,*) word160 ! Line 14      
      READ(7008,*) nele1   ! Line 15
      READ(7008,*) word160 ! Line 16
      READ(7008,*) nele3   ! Line 17
      if (Z_SYMM .eq. 'T') THEN     ! Check for Z_Symmetry
         nele3 = 2.d0 * nele3
         write (9001, *)" Code creates nele3 = 2.d0 * nele3 since Z_SYMM = 'T' "
      endif      
      dx = L1/dble(nele1)  
      write (9001, *)" We set dx = L1/nele1."
      READ(7008,*) word160 ! Line 18
      READ(7008,*) word160 ! Line 19
      READ(7008,*) Csm     ! Line 20     
      READ(7008,*) word160 ! Line 21
      READ(7008,*) Cdm     ! Line 22
      READ(7008,*) word160 ! Line 23
      READ(7008,*) Rom     ! Line 24
      READ(7008,*) word160 ! Line 25
      READ(7008,*) nu      ! Line 26
      READ(7008,*) word160 ! Line 27
      READ(7008,*) csratio ! Line 28
      READ(7008,*) word160 ! Line 29
      READ(7008,*) cdratio ! Line 30
      READ(7008,*) word160 ! Line 31
      READ(7008,*) Roratio ! Line 32
      Mum = Csm*Csm*Rom;   
      cssm = Csm*Csm       
      Cs = Csm/csratio
      Cd = Cdm/cdratio   
      Ro = Rom/Roratio   
      Mu = Cs*Cs*Ro;      
      Muratio = Mum/Mu   
      css = Cs*Cs;  
      factor2 = (1.0d0/(Ro*Cd)+1.0d0/(Rom*Cdm)) ! (c_d/c_s)**2
      factor1 = (Cs/Mu+Csm/Mum)      
      READ(7008,*) word160 ! Line 33
      READ(7008,*) word160 ! Line 34      
      READ(7008,*) Tend    ! Line 35
      READ(7008,*) word160 ! Line 36
      READ(7008,*) betaa   ! Line 37
      dt = betaa*dx/Csm;   ! time step         
      READ(7008,*) word160 ! Line 38
      READ(7008,*) gamma   ! Line 39
      READ(7008,*) word160 ! Line 40
      READ(7008,*) word160 ! Line 41
      READ(7008,*) word160 ! Line 42
      READ(7008,*) Pre_Kernels ! Line 43
      READ(7008,*) word160 ! Line 44
      READ(7008,*) word160 ! Line 45
      READ(7008,*) mus0    ! Line 46
      READ(7008,*) word160 ! Line 47
      READ(7008,*) mur0    ! Line 48
      READ(7008,*) word160 ! Line 49      
      READ(7008,*) DeltaC  ! Line 50
      READ(7008,*) word160 ! Line 51
      READ(7008,*) T0bg    ! Line 52
      READ(7008,*) word160 ! Line 53
      READ(7008,*) word160 ! Line 54      
      READ(7008,*) No_neucles ! Line 55
      If (No_neucles .gt. 0 ) THEN
         Allocate(T0nux11(No_neucles), T0nux31(No_neucles), T0nux12(No_neucles), T0nux32(No_neucles), T0nu(No_neucles,3))
         do i=1, No_neucles
            READ(7008,*) T0nux11(i), T0nux31(i), T0nux12(i), T0nux32(i), T0nu(i,:) ! Line 55 + i
         enddo
      ENDIF
      READ(7008,*) word160 ! Line 56 + No_neucles
      READ(7008,*) word160 ! Line 57 + No_neucles
      READ(7008,*) No_Asp ! Line 58 + No_neucles
      If (No_Asp .gt. 0 ) THEN
         Allocate(Asp_x(No_Asp), Asp_y(No_Asp),Asp_rad(No_Asp),Asp_mus0(No_Asp),Asp_mur0(No_Asp))
         do i=1, No_neucles
            READ(7008,*) Asp_x(i), Asp_y(i), Asp_rad(i), Asp_mus0(i), Asp_mur0(i) ! Line 58 + No_neucles + i
         enddo             
      ENDIF
      READ(7008,*) word160 ! Line 59 + No_neucles + No_Asp
      READ(7008,*) word160 ! Line 60 + No_neucles + No_Asp  
      READ(7008,*) NStations
      IF (NStations .gt. 0) THEN
         Allocate(station_pt(NStations,2)) ! Allocate for output sation point coordinates  
         DO i=1,NStations
            READ(7008,*) station_pt(i,:) ! Line 60 + No_neucles + No_Asp + i
         ENDDO
      ENDIF
   CLOSE(7008)      
   
   
   Nitter = 101
   ntime = ceiling(Tend/dt)
   outinterval1 = ceiling(0.05/dt)
   outinterval2 = ceiling(0.5d0/dt) !12 !6; ! 
   outpos1=  nint(12000.0d0/dx)
   outpos2=  nint(7500.0d0/dx)  ! Position at the interface to compute field quantities
   outpos3=  nint(4500.0d0/dx) 
   Tn0 = T0bg(2) ! For barrier failure stress calculation
   

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Step 1.6: Print input data into a .dat file foir cross-checking 
      WRITE(filename,'(a,a,a)') "./data/",trim(pb_name),"_Kernel_Input.dat"       
      OPEN(unit=7009,file=filename, status='replace', action='write')  
      write(7009,*) 'Printing input perameters for the current job:'
      write(7009,*) 'Solving: ' , trim(TPV_No)
      write(7009,*) 'L1rpt:   ',L1rpt
      write(7009,*) 'L3rpt:   ',L3rpt
      write(7009,*) 'Xr:   ',Xr
      write(7009,*) 'Z_SYMM:   ',Z_SYMM
      write(7009,*) 'L1:   ',L1
      write(7009,*) 'L3:   ',L3        
      write(7009,*) 'nele1:   ', nele1
      write(7009,*) 'nele3:   ', nele3  
      write(7009,*) 'Csm: ', Csm
      write(7009,*) 'Cdm: ', Cdm      
      write(7009,*) 'Rom: ', Rom      
      write(7009,*) 'nu: ', nu
      write(7009,*) 'csratio: ', csratio
      write(7009,*) 'cdratio: ', cdratio
      write(7009,*) 'Roratio: ', Roratio
      write(7009,*) 'dx:   ', dx
      write(7009,*) 'Tend: ', Tend
      write(7009,*) 'betaa:   ', betaa
      write(7009,*) 'gamma:   ', gamma
      write(7009,*) 'dt :   ',dt
      write(7009,*) 'ntime:   ',ntime
      write(7009,*) 'Pre_Kernels  :   ', Pre_Kernels
      write(7009,*) 'mus0  :  ', mus0
      write(7009,*) 'mur0  :   ', mur0
      write(7009,*) 'DeltaC:  ', DeltaC
      write(7009,*) 'T0bg  : ', T0bg  
      write(7009,*) 'No_neucles  : ', No_neucles  
      if ( No_neucles .gt. 0) then
         do i=1, No_neucles
            write (7009,*) T0nux11(i), T0nux31(i), T0nux12(i), T0nux32(i)
            write(7009,*) 'T0nu  : ', T0nu(i,:)
         enddo
      endif      
      
      write(7009,*) 'Tn0   : ', Tn0   
      write(7009,*) 'No_Asp   : ', No_Asp 
      if ( No_Asp .gt. 0) then
         do i=1, No_Asp
            write (7009,*) Asp_x(i), Asp_y(i), Asp_rad(i), Asp_mus0(i), Asp_mur0(i)
         enddo
      endif
      write(7009,*) 'NStations   : ', NStations 
         DO i=1,NStations
            write(7009,*) station_pt(i,:)
         ENDDO   
      close(7009)      
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   ! step 2.1: COMPUTE KERNEL
   n_pts1 = ceiling(dsqrt(2.0d0)*pii*betaa*ntime/(csratio*gamma))
   n_pts3 = ceiling(dsqrt(2.0d0)*pii*betaa*ntime/gamma)
   WRITE(9001,*) 'n_pts1', n_pts1, 'n_pts3', n_pts3, 'gamma',gamma
   flush(9001)
   allocate(m11u(0:n_pts1),tu(0:n_pts1),m11d(0:n_pts3),td(0:n_pts3),m22u(0:n_pts1),m22d(0:n_pts3),m12u(0:n_pts1),m12d(0:n_pts3))   
 
 
   
 ! Printing the information of Pre-computed Kernels 
   Write (9001,*)'Pre-Kernels Set to:' , Pre_Kernels
   Write (9001,*)'Meaning'
   flush(9001)   
   if (Pre_Kernels == 1) then
      Write (9001,*)'Pre-computed Kernels set as available.'
      Write(9001,*)'Make sure Kernels are pre-computed and stored in ./',trim(pb_name),'_Kernels directory.'
      flush(9001)   
   endif
   if (Pre_Kernels == 0) then
      Write (9001,*)'Pre-computed Kernels set as not available.'
      Write (9001,*)'Computing the kernels'
      flush(9001)   
   endif 
 
   ! This loop calculates or retrives the kernels as per the input value

   ! This loop calculates or retrives the kernels as per the input value   
   if (Pre_Kernels == 1) then
      Write (9001,*)'Loading u kernel' 
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m11u.dat"
      OPEN(91,file=trim(filename))
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m11d.dat"       
      OPEN(92,file=trim(filename))      
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m22u.dat"       
      OPEN(93,file=trim(filename))      
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m22d.dat"       
      OPEN(94,file=trim(filename))      
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m12u.dat"       
      OPEN(95,file=trim(filename))      
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m12d.dat"       
      OPEN(96,file=trim(filename))      
      read(91,*)(m11u(i),i=1,n_pts1)
      Write (9001,*) 'm11u read', size(m11u), m11u(1)
      read(92,*)(m11d(i),i=1,n_pts3)
      Write (9001,*) 'm11d read', size(m11d), m11d(1)
      read(93,*)(m22u(i),i=1,n_pts1)
      Write (9001,*) 'm22u read', size(m22u), m22u(1)
      read(94,*)(m22d(i),i=1,n_pts3)
      Write (9001,*) 'm22d read', size(m22d), m22d(1)
      read(95,*)(m12u(i),i=1,n_pts1)
      Write (9001,*) 'm12u read', size(m12u), m12u(1)
      read(96,*)(m12d(i),i=1,n_pts3)
      Write (9001,*) 'm12d read', size(m12d), m12d(1)
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      close(96)
   elseif (Pre_Kernels == 0) then
      CALL kern11(n_pts1,gamma,m11u,tu)
      ! open(1001,file='./Kernels/Kernel_m11u.dat')      
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m11u.dat"       
      OPEN(1001,file=trim(filename))      
      do i = 1,n_pts1
         write(1001,*) m11u(i)
      enddo
      Write (9001,*) 'm11u read', size(m11u), m11u(1)
      close(1001)
      
      CALL kern11(n_pts3,gamma,m11d,td)
      ! open(2001,file='./Kernels/Kernel_m11d.dat')
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m11d.dat"       
      OPEN(2001,file=trim(filename))      
      do i = 1,n_pts3
         write(2001,*) m11d(i)
      enddo
      Write (9001,*) 'm11d read', size(m11d), m11d(1)
      close(2001)
      
      CALL kern22(n_pts1,gamma,m22u,tu)
      ! open(1002,file='./Kernels/Kernel_m22u.dat')
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m22u.dat"       
      OPEN(1002,file=trim(filename))           
      do i = 1,n_pts1
         write(1002,*) m22u(i)
      enddo
      Write (9001,*) 'm22u read', size(m22u), m22u(1)
      close(1002)
      
      CALL kern22(n_pts3,gamma,m22d,td)
      ! open(2002,file='./Kernels/Kernel_m22d.dat')
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m22d.dat"       
      OPEN(2002,file=trim(filename))      
      do i = 1,n_pts3
         write(2002,*) m22d(i)
      enddo
      Write (9001,*) 'm22d read', size(m22d), m22d(1)
      close(2002)
      
      CALL kern12(n_pts1,gamma,m12u,tuu)
      ! open(1012,file='./Kernels/Kernel_m12u.dat')
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m12u.dat"       
      OPEN(1012,file=trim(filename))          
      do i = 1,n_pts1
         write(1012,*) m12u(i) 
      enddo
      Write (9001,*) 'm12u read', size(m12u), m12u(1)
      close(1012)
      
      CALL kern12(n_pts3,gamma,m12d,tdd)
      ! open(2012,file='./Kernels/Kernel_m12d.dat')
      WRITE(filename,'(a,a,a,a)') "./",trim(pb_name),'_Kernels',"/Kernel_m12d.dat"       
      OPEN(2012,file=trim(filename))            
      do i = 1,n_pts3
         write(2012,*) m12d(i)
      enddo
      Write (9001,*) 'm12d read', size(m12d), m12d(1)
      close(2012)
   else
      Write (9001,*)'Input Error'
      Write (9002,*)'Neither Pre_Kernels == 0 Nor Pre_Kernels == 1'
      stop
   endif
   m22u(0) = 0.0d0
   m22d(0) = 0.0d0
   m11u(0) = 0.0d0
   m11d(0) = 0.0d0
   flush(9001)
   flush(9002)   
   

   Write (9001,*) '####### KERNEL CHECK TEST DONE ########'
   ! GET THE END TIME
   CALL CPU_TIME(END_TIME)
   ! CALCULATE THE ELAPSED TIME
   ELAPSED_TIME = END_TIME - START_TIME
   ! OUTPUT THE ELAPSED TIME
   Write (9001,*) "ELAPSED TIME:", ELAPSED_TIME, "SECONDS"
   close(9001)
   close(9002)
END PROGRAM Kernel_Check