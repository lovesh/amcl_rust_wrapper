# amcl_rust_wrapper

Wraps some parts of [AMCL](https://github.com/miracl/amcl) to provide a nice abstraction to work with finite field elements and group elements when working with elliptic curves.   
Overloads +, -, *, +=, -= to use with field as well as group elements.  
Provides abstraction for creating vectors of field elements or group elements and then scale, add, subtract, take inner product or Hadamard product.   
Additionally, implements some extra algoirthms like variable time scalar multiplication using wNAF, constant time and variable time multi-scalar multiplication, batch (simultaneous) inversion and Barrett reduction.
