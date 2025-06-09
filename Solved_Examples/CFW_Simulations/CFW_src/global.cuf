! module GLOBAL starting from here

MODULE global

  REAL*8 :: T, pi, limita, limitb
  REAL*8 :: yy1, yy2, yyr, yr, aa, a, fyr, gyr
CONTAINS

  SUBROUTINE set_constants()
  IMPLICIT NONE 
  pi  = 4.0d0*datan(1.0d0)
  yyr = 2.0d0-(2.0d0/dsqrt(3.0d0))
  yy1 = 2.0d0+(2.0d0/dsqrt(3.0d0))
  yy2 = 4.0d0
  yr  = dsqrt(yyr)
  aa  = 3.0d0
  a   = dsqrt(aa)
  fyr = f(yr)
  gyr = g(yr)
  END SUBROUTINE set_constants

  FUNCTION f(y)
  IMPLICIT NONE
  REAL*8 :: f, y, yy, neum, denom
  yy = y*y
  neum  = (4.0d0*(1.0d0-yy)*(dsqrt(1.0d0-(yy/aa))))+((dsqrt(1.0d0-yy))*(2.0d0-yy)**2)
  denom = (y+yr)*(yy-yy1)*(yy-yy2)
  f     = (2.0d0/pi)*(neum/denom)
  END FUNCTION f

  FUNCTION g1(y)
  IMPLICIT NONE
  REAL*8 :: g1, y
  g1 = -1*((f(y)-fyr)/(y-yr))*y*dsin(y*T)
  END FUNCTION g1

  FUNCTION g2(y)
  IMPLICIT NONE
  REAL*8 :: g2, y, yy, neum, denom
  yy    = y*y
  neum  = (4.0d0*((1-yy)))*dsqrt(1.0d0-yy/aa)
  denom = (yy-yyr)*(yy-yy1)*(yy-yy2)
  g2    = (-2.0d0/pi)*(neum/denom)*y*dsin(y*T)
  END FUNCTION g2
  
  !define function cosx


  FUNCTION g(y)
  IMPLICIT NONE
  REAL*8 :: g, y, yy, neum, denom
  yy = y*y
  neum  = (4.0d0*(dsqrt(1.0d0-yy))*(1.0d0-(yy/aa)))+((dsqrt(1.0d0-(yy/aa)))*(2.0d0-yy)**2)
  denom = (y+yr)*(yy-yy1)*(yy-yy2)
  g     = (2.0d0/pi)*(neum/denom)
  END FUNCTION g


  FUNCTION f1(y)
  IMPLICIT NONE
  REAL*8 :: f1, y
  f1 = -1*((g(y)-gyr)/(y-yr))*y*dsin(y*T)
  END FUNCTION f1


  FUNCTION f2(y)
  IMPLICIT NONE
  REAL*8 :: f2, y, yy, neum, denom
  yy    = y*y
  neum  = ((2.0d0-yy)**2)*dsqrt(1.0d0-yy/aa)
  denom = (yy-yyr)*(yy-yy1)*(yy-yy2)
  f2    = (-2.0d0/pi)*(neum/denom)*y*dsin(y*T)
  END FUNCTION f2  

 
 FUNCTION b2(y)
  IMPLICIT NONE
  REAL*8:: b2, y, neum, denom, yy
  yy = y*y
  neum = (2.0d0*y*(2.0d0-yy)*dsqrt(yy-1.0d0)*dsqrt(1.0d0-(yy/aa)))
  denom = (yyr-yy)*(yy1-yy)*(yy2-yy)
  b2 = (neum/denom)*dcos(y*T)
 END FUNCTION b2
 

FUNCTION cosx(x)
 IMPLICIT NONE
  REAL*8:: cosx
  REAL*8, INTENT(IN)::x
 cosx=(dcos(x))/x
END FUNCTION cosx

!define function sinx


FUNCTION sinx(x)
 IMPLICIT NONE
  REAL*8:: sinx
  REAL*8, INTENT(IN)::x
 sinx=(dsin(x))/x
END FUNCTION sinx

END MODULE global