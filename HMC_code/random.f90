function r_normal ()
!	random number generator
!	sampling from the normal distribution
!	unform random number obtained from 

	implicit none

	real(kind=8)::u1,u2
	real(kind=8)r_normal
	real(kind=8),parameter::pi = 3.141592653589793d0

	call random_number(u1)
	call random_number(u2)
!	Box-Muller
	r_normal=sqrt(-2.d0*log(u1))*cos(2.d0*pi*u2)

	return
end
