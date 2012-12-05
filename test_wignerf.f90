program test_wigner

	use wignercoupling

	implicit none

	integer, parameter :: dp = kind(1.0d0)
	real(dp), dimension(:), allocatable :: temp
	integer :: i, j

	write (*,*) wigner3jf(2,2,0,0,0,0), -sqrt(1.0_dp/3.0_dp)
	write (*,*) wigner3jf(8,10,6,4,-4,0), sqrt(1.0_dp/2145.0_dp)
	write (*,*) wigner3jf(6,6,6,2,-2,0), -sqrt(1.0_dp/42.0_dp)
	write (*,*) wigner3jf(3,9,6,1,3,-4), 0.5_dp*sqrt(3.0_dp/35.0_dp)

	temp =  wigner3jvectf(2,2,0,0)
	write (*,*) temp
	write (*,*) -sqrt(1.0_dp/3.0_dp), 0.0_dp, sqrt(2.0_dp/15.0_dp)

	write (*,*) wigner3jvectf(2*2,2*2,(-2+1)*2,2*2)
	write (*,*) wigner3jvectf(2*2,2*2,(2-1)*2,-2*2)
	write (*,*) wigner3jvectf(200*2,2*2,0,0)

end program test_wigner
