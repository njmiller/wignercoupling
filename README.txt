Compiles a shared (static) library that calculates Wigner 3j/6j/9j symbols

This calculation is based off of a recursive formula and not factorials so it should be able to calculate the Wigner values for large numbers without overflows.

Bindings for IDL and Fortran were also written.

This was written in the mid 2000s so there may be more useful libraries now. I am not sure how much this version was tested as I rewrote it in Python because I mainly use that for code now and tested that version much more.
