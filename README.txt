Compiles a shared (static) library that calculates Wigner 3j/6j/9j symbols

Can't really use OpenMP because what goes on in loop depends on what previously 
went on

Things to do:
	Convert everything to C
	Test and fix Wigner 6j
	Test and fix Wigner 9j
	Write a DLM so everything can be called from IDL (why I am converting everything to C).
