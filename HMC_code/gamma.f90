!---------------
FUNCTION random_gamma(s, b, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     N.B. This version is in `double precision' and includes scaling

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).
!     B = Scale parameter

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: s, b
LOGICAL, INTENT(IN)    :: first
REAL (dp)              :: fn_val

! Local parameters
REAL (dp), PARAMETER  :: one = 1.0_dp, zero = 0.0_dp

IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s >= one) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < one) THEN
  fn_val = random_gamma2(s, first)
END IF

! Now scale the random variable
fn_val = b * fn_val
RETURN

CONTAINS


FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO GAMMA**(S-1)*EXP(-GAMMA),
! BASED UPON BEST'S T DISTRIBUTION METHOD

!     S = SHAPE PARAMETER OF DISTRIBUTION
!          (1.0 < REAL)

REAL (dp), INTENT(IN)  :: s
LOGICAL, INTENT(IN)    :: first
REAL (dp)              :: fn_val

!     Local variables
REAL (dp)             :: d, r, g, f, x
REAL (dp), SAVE       :: b, h
REAL (dp), PARAMETER  :: sixty4 = 64.0_dp, three = 3.0_dp, pt75 = 0.75_dp,  &
                         two = 2.0_dp, half = 0.5_dp

IF (s <= one) THEN
  WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  b = s - one
  h = SQRT(three*s - pt75)
END IF

DO
  CALL RANDOM_NUMBER(r)
  g = r - r*r
  IF (g <= zero) CYCLE
  f = (r - half)*h/SQRT(g)
  x = b + f
  IF (x <= zero) CYCLE
  CALL RANDOM_NUMBER(r)
  d = sixty4*g*(r*g)**2
  IF (d <= zero) EXIT
  IF (d*x < x - two*f*f) EXIT
  IF (LOG(d) < two*(b*LOG(x/b) - f)) EXIT
END DO
fn_val = x

RETURN
END FUNCTION random_gamma1



FUNCTION random_gamma2(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

REAL (dp), INTENT(IN)  :: s
LOGICAL, INTENT(IN)    :: first
REAL (dp)              :: fn_val

!     Local variables
REAL (dp)             :: r, x, w
REAL (dp), SAVE       :: a, p, c, uf, vr, d
REAL (dp), PARAMETER  :: vsmall = EPSILON(one)

IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2

END FUNCTION random_gamma
