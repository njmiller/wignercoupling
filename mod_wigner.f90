module wignercoupling
!Module that is an interface to the C routines

	use iso_c_binding

	implicit none

	integer, parameter :: sp = kind(1.0)
	integer, parameter :: dp = kind(1.0d0)

contains

	pure function wigner3jf(tj1, tj2, tj3, tm1, tm2, tm3)
	
		real(dp) :: wigner3jf
		integer(sp), intent(in) :: tj1, tj2, tj3, tm1, tm2, tm3
		integer(c_int) :: tj1p, tj2p, tj3p, tm1p, tm2p, tm3p

		interface
			real(kind=c_double) function Wigner3j(tj1,tj2,tj3,tm1,tm2,tm3) bind(c)
				use iso_c_binding
				integer(kind=c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
			end function Wigner3j
		end interface

		tj1p = tj1
		tj2p = tj2
		tj3p = tj3
		tm1p = tm1
		tm2p = tm2
		tm3p = tm3
		wigner3jf = Wigner3j(tj1p,tj2p,tj3p,tm1p,tm2p,tm3p)

	end function wigner3jf

end module wignercoupling
