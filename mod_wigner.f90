module wignercoupling
!Module that is an interface to the C routines

	use iso_c_binding

	implicit none

	private

	integer, parameter :: sp = kind(1.0)
	integer, parameter :: dp = kind(1.0d0)

	public :: wigner3jf, clebschgordonf, racahvf, wigner6jf

contains
	
	pure function wigner3jf(tj1, tj2, tj3, tm1, tm2, tm3)
	
		real(dp) :: wigner3jf
		integer(sp), intent(in) :: tj1, tj2, tj3, tm1, tm2, tm3
		integer(c_int) :: tj1p, tj2p, tj3p, tm1p, tm2p, tm3p

		interface
			pure function wigner3j(tj1,tj2,tj3,tm1,tm2,tm3) bind(c)
				use iso_c_binding
				real(kind=c_double) :: wigner3j
				integer(kind=c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
			end function wigner3j
		end interface

		tj1p = tj1
		tj2p = tj2
		tj3p = tj3
		tm1p = tm1
		tm2p = tm2
		tm3p = tm3
		wigner3jf = wigner3j(tj1p,tj2p,tj3p,tm1p,tm2p,tm3p)

	end function wigner3jf

	pure function clebschgordonf(tj1,tj2,tj3,tm1,tm2,tm3)
		!Even though it is trivial to calculate this from Wigner 3j symbols
		!we will call the C function so that it is only implemented in one place
		real(dp) :: clebschgordonf
		integer(sp), intent(in) :: tj1, tj2, tj3, tm1, tm2, tm3
		integer(c_int) :: tj1p, tj2p, tj3p, tm1p, tm2p, tm3p

		interface
			pure function clebschgordon(tj1,tj2,tj3,tm1,tm2,tm3) bind(c)
				use iso_c_binding
				real(kind=c_double) :: clebschgordon
				integer(kind=c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
			end function clebschgordon
		end interface

		tj1p = tj1
		tj2p = tj2
		tj3p = tj3
		tm1p = tm1
		tm2p = tm2
		tm3p = tm3
		clebschgordonf = clebschgordon(tj1p,tj2p,tj3p,tm1p,tm2p,tm3p)

	end function clebschgordonf

	pure function racahvf(tj1,tj2,tj3,tm1,tm2,tm3)	
		!Even though it is trivial to calculate this from Wigner 3j symbols
		!we will call the C function so that it is only implemented in one place
		real(dp) :: racahvf
		integer(sp), intent(in) :: tj1, tj2, tj3, tm1, tm2, tm3
		integer(c_int) :: tj1p, tj2p, tj3p, tm1p, tm2p, tm3p

		interface
			pure function racahv(tj1,tj2,tj3,tm1,tm2,tm3) bind(c)
				use iso_c_binding
				real(kind=c_double) :: racahv
				integer(kind=c_int), value :: tj1, tj2, tj3, tm1, tm2, tm3
			end function racahv
		end interface

		tj1p = tj1
		tj2p = tj2
		tj3p = tj3
		tm1p = tm1
		tm2p = tm2
		tm3p = tm3
		racahvf = racahv(tj1p,tj2p,tj3p,tm1p,tm2p,tm3p)

	end function racahvf
	
	pure function wigner6jf(tj1,tj2,tj3,tk1,tk2,tk3)

		real(dp) :: wigner6jf
		integer(sp), intent(in) :: tj1, tj2, tj3, tk1, tk2, tk3
		integer(c_int) :: tj1p, tj2p, tj3p, tk1p, tk2p, tk3p

		interface
			pure function wigner6j(tj1,tj2,tj3,tk1,tk2,tk3) bind(c)
				use iso_c_binding
				real(kind=c_double) :: wigner6j
				integer(kind=c_int), value :: tj1, tj2, tj3, tk1, tk2, tk3
			end function wigner6j
		end interface

		tj1p = tj1
		tj2p = tj2
		tj3p = tj3
		tk1p = tk1
		tk2p = tk2
		tk3p = tk3

		wigner6jf = wigner6j(tj1p,tj2p,tj3p,tk1p,tk2p,tj3p)

	end function wigner6jf

end module wignercoupling
