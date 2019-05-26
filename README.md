# amcl_rust_wrapper

Wraps some parts of [AMCL](https://github.com/miracl/amcl) to provide a nice abstraction to work with finite field elements and group elements when working with elliptic curves.   
Overloads +, -, *, +=, -= to use with field as well as group elements.  
Additionaly, implements some extra algoirthms like variable time scalar multiplication using wNAF, constant time and variable time multi-scalar multiplication, batch (simultaneous) inversion and Barrett reduction.
