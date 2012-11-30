Compiles a shared (static) library that calculates Wigner 3j/6j/9j symbols

The Wigner 3j functions seem to be working. I tested all the different cases and 
got the same result as an online calculator that gave me symbolic results.

Can't really use OpenMP because what goes on in loop depends on the previous 
iteration.

Everything is written in C so that I can easily create bindings to other languages. I 
currently have written a binding to Fortran.

Things to do:
	Test and fix Wigner 6j
	Test and fix Wigner 9j
	Write a DLM so everything can be called from IDL (why I am converting everything to C).
	Write Cython module
	Write D module (it may not be needed since all that is needed is a function prototype)
