
! module kernel starting from here

MODULE kernel 

 CONTAINS

  !inserting subroutine for trapezoidal
  SUBROUTINE integral(func,n,ofunc,x0,xn)
    IMPLICIT NONE
    REAL*8, EXTERNAL::func
    INTEGER, INTENT(IN)::n
    REAL*8, INTENT(IN)::x0,xn
    REAL*8, INTENT(OUT)::ofunc
    INTEGER ::i
    REAL*8 ::h,sum
    h=(xn-x0)/n
    sum=(func(x0)+func(xn))/2.0d0
! The OpenMP is utilized to execute the loop, enabling it to operate concurrently on multiple CPUs.
!$OMP PARALLEL DO REDUCTION(+:sum)
    DO i=1, n-1
    sum=sum+func(x0+(i*h))
    END DO
!$OMP END PARALLEL DO
    ofunc=sum*h
END SUBROUTINE integral 


! This subroutine calculate kernel, m11(k,t), in convolution term.
  SUBROUTINE kern11(npt,gamma,m11,tt)
  USE global 
  
  IMPLICIT NONE
  INTEGER :: start, fini, i, iweigh, npt
  REAL*8  :: gamma
  REAL*8, DIMENSION(:) :: m11,tt
  REAL*8, ALLOCATABLE, DIMENSION(:) :: t1, t2, t3, t4, t5, lt, t01, t02, t11, t12
  CALL set_constants()
  allocate(t1(npt),t01(npt),t02(npt),t2(npt),t11(npt),t12(npt),t3(npt),t4(npt),t5(npt),lt(npt))
   
   do i=1,npt
    tt(i)=dfloat(i)*gamma
    T = tt(i)
    
    !done repair of f1 to split up
    call integral(g1, 100000, t01(i), 0.0d0, 0.919400686d0)
    call integral(g1, 100000, t02(i), 0.919402686d0, 1.0d0)
    !done for func f2
    call integral(g2, 100000, t2(i), 1.0d0, a)
    
    t1(i)=t01(i)+t02(i)
   
    !space for cos integration limits with 0.000000001d0 onwards for sin and cos all
    limita=(1.0d0-yr)*T
    limitb=yr*T
    
    call integral(cosx, 100000, t11(i), limita, limitb)
    call integral(sinx, 100000, t12(i), -limitb, limita)
   
    t3(i)=fyr*yr*dsin(yr*T)*t11(i) 
    t4(i)=-fyr*yr*dcos(yr*T)*t12(i)
    if (dabs(T) .le. 1e-6) then
      t5(i) = 0.0d0
    else
      t5(i) =  fyr*(dcos(T)-1.d0)/T
    endif
    lt(i)  = -fyr*yr*pi*dcos(yr*T)
    m11(i) =  t1(i) + t2(i) + t3(i) + t4(i) + t5(i)
  enddo 

deallocate (t1, t2, t3, t4, t5, lt, t11, t12)

  END SUBROUTINE kern11

! This subroutine calculate kernel, m22(k,t), in convolution term.
  SUBROUTINE kern22(npt,gamma,m22,tt)
  USE global
  
  IMPLICIT NONE
  INTEGER :: start, fini, i, iweigh, npt
  REAL*8  :: gamma
  REAL*8, DIMENSION(:) :: m22,tt
  REAL*8, ALLOCATABLE, DIMENSION(:) :: t1, t2, t3, t4, t5, lt, t01, t02, t11, t12
  CALL set_constants()
  allocate(t1(npt),t01(npt),t02(npt),t2(npt),t11(npt),t12(npt),t3(npt),t4(npt),t5(npt),lt(npt))

  do i=1,npt
    tt(i)=dfloat(i)*gamma
    T    = tt(i)
    
    !done repair of f1 to split up
    call integral(f1, 100000, t01(i), 0.0d0, 0.919400686d0)
    call integral(f1, 100000, t02(i), 0.919402686d0, 1.0d0)
    !done for func f2
    call integral(f2, 100000, t2(i), 1.0d0, a)
    
    t1(i)=t01(i)+t02(i)
   
    !space for cos integration limits with 0.000000001d0 onwards for sin and cos all
    limita=(1.0d0-yr)*T
    limitb=yr*T
    
    call integral(cosx, 100000, t11(i), limita, limitb)
    call integral(sinx, 100000, t12(i), -limitb, limita)
   
    t3(i)=gyr*yr*dsin(yr*T)*t11(i) 
    t4(i)=-gyr*yr*dcos(yr*T)*t12(i)
    if (dabs(T) .le. 1e-6) then
      t5(i) = 0.0d0
    else
      t5(i) =  gyr*(dcos(T)-1.d0)/T
    endif
    lt(i)  = -gyr*yr*pi*dcos(yr*T)
    m22(i) =  t1(i) + t2(i) + t3(i) + t4(i) + t5(i)
  enddo

deallocate (t1, t2, t3, t4, t5, lt,  t11, t12)
  END SUBROUTINE kern22

! This subroutine calculate kernel, m12(k,t), in convolution term.
SUBROUTINE kern12(npt,gamma,m12,tt)
  USE global
  
  IMPLICIT NONE
  INTEGER :: i, npt
 
  REAL*8  :: gamma, const, neum, neum1, neum2, denom
  REAL*8, ALLOCATABLE, DIMENSION(:) :: m12,tt
  REAL*8, ALLOCATABLE, DIMENSION(:) :: t1, t2, t3
 
  CALL set_constants()
  allocate(t1(0:npt),t2(0:npt),t3(0:npt))
  allocate(m12(0:npt),tt(0:npt))
  !constant term of second term done here
  
  neum1 = (-8.0d0*yyr/aa)+(8.0d0/aa)+(-1.0d0*yyr*yyr)+(-4.0d0)+(6.0d0*yyr)
  neum2 = (2.0d0*(2.0d0-yyr)*dsqrt(1.0d0-yyr)*dsqrt(1.0d0-(yyr/aa)))
  neum  = neum1+neum2
  denom = (yy1-yyr)*(yy2-yyr) 
  const = neum/denom
  
 do i=0,npt
   tt(i) = dfloat(i)*gamma
    T    = tt(i)
   t1(i) = const*dcos(yr*T)
   
   call integral(b2, 100000, t2(i), 1.0d0, a)
   
    t3(i) = (-2.0d0/pi)*t2(i)
   m12(i) = -(t1(i) + t3(i))     
 
 enddo 
 END SUBROUTINE kern12
  
END MODULE kernel