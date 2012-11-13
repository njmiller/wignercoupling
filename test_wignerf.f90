program test_wigner

	use wignercoupling

	implicit none

	integer, parameter :: dp = kind(1.0d0)

	real(dp) :: wigval

	write (*,*) wigner3jf(2,2,0,0,0,0), -0.577350
	write (*,*) wigner3jf(8,10,6,4,-4,0), 0.0215917
	write (*,*) wigner3jf(6,6,6,2,-2,0), -0.1543033
	write (*,*) wigner3jf(3,9,6,1,3,-4), 0.146385

end program test_wigner
