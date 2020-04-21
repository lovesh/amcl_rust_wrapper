# amcl_rust_wrapper

- Wraps some parts of [AMCL](https://github.com/miracl/amcl) to provide a nice abstraction to work with finite field elements and group elements when working with elliptic curves.   
- Overloads +, -, *, +=, -= to use with field as well as group elements.  The overloaded operators correspond to constant time methods. But for scalar multiplication, variable time algorithms are present but can be used by calling methods only. 
- Provides abstraction for creating vectors of field elements or group (elliptic curve points) elements and then scale, add, subtract, take inner product or Hadamard product.
- Supports creating univariate polynomials of field elements and doing arithmetic on them.
- Some of the operations on vectors and polynomials are parallelized using [rayon](https://github.com/rayon-rs/rayon).
- Serialization support using [serde](https://github.com/serde-rs/json).
- Field and group elements are cleared when dropped. Using [zeroize](https://crates.io/crates/zeroize).     
- Additionally, implements some extra algorithms like variable time scalar multiplication using wNAF, constant time and variable time multi-scalar multiplication, batch (simultaneous) inversion and Barrett reduction.

## Building
The wrapper has to be built by enabling any one of the mentioned curve as a feature.   
To build for BLS-381 curve, use 

```
cargo build --no-default-features --features bls381
```

Similarly, to build for secp256k1, use
```
cargo build --no-default-features --features secp256k1
```

To run tests for secp256k1, use
```
cargo test --no-default-features --features secp256k1
```

To use it as dependency crate, use the name of the curve as a feature. eg. for BLS12-381 curve, use
```
[dependencies.amcl_wrapper]
version = "0.3.4"
default-features = false
features = ["bls381"]
```

Note that only one curve can be used at a time so the code only works with one feature.

## Benchmarking
There are tests for various operations which print the time taken to do those ops. They are prefixed with `timing`*[]: 
To run them use

```
cargo test --release --no-default-features --features <curve name> -- --nocapture timing
```

## Examples
1. Create some random field elements or group elements and do some basic additions/subtraction/multiplication 
```rust
let a = FieldElement::random();
let b = FieldElement::random();
let neg_b = -b;
assert_ne!(b, neg_b);
let neg_neg_b = -neg_b;
assert_eq!(b, neg_neg_b);
assert_eq!(b+neg_b, FieldElement::zero());
        
let c = a + b;
let d = a - b;
let e = a * b;
let mut sum = FieldElement::zero();
sum += c;
sum += d;
```

```rust
// Compute inverse of element
let n = a.invert();

// Faster but constant time computation of inverse of several elements. `inverses` will be a vector of inverses of elements and `pr` will be the product of all inverses.
let (inverses, pr) = FieldElement::batch_invert(vec![a, b, c, d].as_slice());  
```

```rust
// Compute a width 5 NAF.
let a_wnaf = a.to_wnaf(5);
```

```rust
// G1 is the elliptic curve sub-group over the prime field
let x = G1::random();
let y = G1::random();
let neg_y = -y;
assert_ne!(y, neg_y);
let neg_neg_y = -neg_y;
assert_eq!(y, neg_neg_y);
assert_eq!(y+neg_y, G1::identity());

let z = x + y;
let z1 = x - y;
let mut sum_1 = G1::identity();
sum_1 += z;
sum_1 += z1;
```

```rust
// G2 is the elliptic curve sub-group over the prime extension field
let x = G2::random();
let y = G2::random();
let neg_y = -y;
assert_ne!(y, neg_y);
let neg_neg_y = -neg_y;
assert_eq!(y, neg_neg_y);
assert_eq!(y+neg_y, G2::identity());

let z = x + y;
let z1 = x - y;
let mut sum_1 = G2::identity();
sum_1 += z;
sum_1 += z1;

// To check that G1 or G2 have correct order, i.e. the curve order
let x = G1::random();
assert!(x.has_correct_order());
let y = G2::random();
assert!(y.has_correct_order());
```

Mutating versions of the above operations like addition/subtraction/negation/inversion are present but have to be called as methods like `b.negate()`

2. Scalar multiplication
```rust
let a = FieldElement::random();
let g = G1::generator();  // the group's generator
// constant time scalar multiplication
let m = g * a;
// variable time scalar multiplication using wNAF
let n = g.scalar_mul_variable_time(&a);
```

3. Map an arbitrary size message to a field element or group element by hashing the message internally.
```
let msg = "Some message";
let a = FieldElement::from_msg_hash(msg.as_bytes());
let b = G1::from_msg_hash(msg.as_bytes());
```

4. Create vectors of field elements and do some operations
```rust
// creates a vector of size 10 with all elements as 0
let mut a = CurveOrderElementVector::new(10);
// Add 2 more elements to the above vector
a.push(FieldElement::random());
a.push(FieldElement::random());

a[0];       // 0th element of above vector
a[1];       // 1st element of above vector

a.len();    // length of vector
a.sum();    // sum of elements of vector 
```

```rust
// Return a Vandermonde vector of a given field element, i.e. given element `k` and size `n`, return vector as `vec![1, k, k^2, k^3, ... k^n-1]`
let k = FieldElement::random();  
let van_vec = CurveOrderElementVector::new_vandermonde_vector(&k, 5);
```

```rust
// creates a vector of size 10 with randomly generated field elements
let rands: Vec<_> = (0..10).map(|_| FieldElement::random()).collect();

// an alternative way of creating vector of size 10 of random field elements
let rands_1 = CurveOrderElementVector::random(10);
```

```rust
// Compute new vector as sum of 2 vectors. Requires vectors to be of equal length. 
let sum_vec = rands.plus(&rands_1);
// Compute new vector as difference of 2 vectors. Requires vectors to be of equal length.
let diff_vec = rands.minus(&rands_1);
// Return the scaled vector of the given vector by a field element `n`, i.e. multiply each element of the vector by that field element
let n = FieldElement::random();
let scaled_vec = rands.scaled_by(&n);
// Scale the vector itself
let mut rands_2 = rands_1.clone();
rands_2.scale(&n);
```

```rust
// Compute inner product of 2 vectors. Requires vectors to be of equal length.
let ip = rands.inner_product(&rands_1);

// Compute Hadamard product of 2 vectors. Requires vectors to be of equal length.
let hp = rands.hadamard_product(&rands_1);
```

5. Create vectors of group elements and do some operations
```rust
// creates a vector of size 10 with all elements as 0
let mut a = G1Vector::new(10);
// Add 2 more elements to the above vector
a.push(G1::random());
a.push(G1::random());


// creating vector of size 10 of random group elements
let rands: Vec<_> = (0..10).map(|_| G1::random()).collect();
let rands_1 = G1Vector::random(10);
// Compute new vector as sum of 2 vectors. Requires vectors to be of equal length. 
let sum_vec = rands.plus(&rands_1);
// Compute new vector as difference of 2 vectors. Requires vectors to be of equal length.
let diff_vec = rands.minus(&rands_1);
```

```rust
// Compute inner product of a vector of group elements with a vector of field elements.    
// eg. given a vector of group elements and field elements, G and F respectively, compute G[0]*F[0] + G[1]*F[1] + G[2]*F[2] + .. G[n-1]*F[n-1]   
// requires vectors to be of same length
let g = G1Vector::random(10);
let f = CurveOrderElementVector::random(10);

// Uses constant time multi-scalar multiplication `multi_scalar_mul_const_time` underneath. 
let ip = g.inner_product_const_time(&f);

// Uses variable time multi-scalar multiplication `multi_scalar_mul_var_time` underneath.
let ip1 = g.inner_product_var_time(&f);

// If lookup tables are already constructed, `multi_scalar_mul_const_time_with_precomputation_done` and `multi_scalar_mul_var_time_with_precomputation_done` can be used for constant and variable time multi-scalar multiplication
```

5. Pairing support. Ate pairing is supported with target group `GT`  
```rust
let g1 = G1::random();
let g2 = G2::random();
// compute reduced ate pairing for 2 elements, i.e. e(g1, g2)
let gt = GT::ate_pairing(&g1, &g2);

// multiply target group elements
let h1 = G1::random();
let h2 = G2::random();
let ht = GT::ate_pairing(&h1, &h2);
let m = GT::mul(&gt, &ht);

// compute reduced ate pairing for 4 elements, i.e. e(g1, g2) * e (h1, h2)
let p = GT::ate_2_pairing(&g1, &g2, &h1, &h2);

// compute reduced ate multi-pairing. Takes a vector of tuples of group elements G1 and G2 as Vec<(&G1, &G2)>
let e = GT::ate_multi_pairing(tuple_vec);

// Raise target group element to field element (GT^f)
let r = FieldElement::random();
let p = gt.pow(&r);
```

6. Serialization
```rust
// Convert to and from hex.
let hex_repr = a.to_hex();
let a_recon = FieldElement::from_hex(hex_repr).unwrap();    // Constant time conversion

// Convert to and from hex.
let hex_repr = x.to_hex();
let x_recon = G1::from_hex(hex_repr).unwrap();    // Constant time conversion

// Convert to and from hex.
let hex_repr = y.to_hex();
let y_recon = G2::from_hex(hex_repr).unwrap();    // Constant time conversion
```

7. Univariate polynomials
```rust
// Create a univariate zero polynomial of degree `d`, i.e. the polynomial will be 0 + 0*x + 0*x^2 + 0*x^3 + ... + 0*x^d 
let poly = UnivarPolynomial::new(d);
assert!(poly.is_zero());

// Create a polynomial from field elements as coefficients, the following polynomial will be c_0 + c_1*x + c_2*x^2 + c_3*x^3 + ... + c_d*x^d
let coeffs: Vec<FieldElement> = vec![c_0, c_1, ... coefficients for smaller to higher degrees ..., c_d];
let poly1 = UnivarPolynomial(CurveOrderElementVector::from(coeffs));

// Create a polynomial of degree `d` with random coefficients 
let poly2 = UnivarPolynomial::random(d);

// Create a polynomial from its roots  
let poly3 = UnivarPolynomial::new_with_roots(roots);

// Create a polynomial by passing the coefficients to a macro
let poly4 = univar_polynomial!(
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::from(87u64),
            -FieldElement::one(),
            FieldElement::from(300u64)
        );

// A polynomial can be evaluated at a field element `v` 
let res: FieldElement = poly1.eval(v);

// Polynomials can be added, subtracted, multiplied or divided to give new polynomials
let sum = UnivarPolynomial::sum(&poly1, &poly2);
// Or use operator overloading
let sum = &poly1 + &poly2;
let diff = UnivarPolynomial::difference(&poly1, &poly2);
// Or use operator overloading
let diff = &poly1 - &poly2;
let product = UnivarPolynomial::multiply(&poly1, &poly2);
// Or use operator overloading
let product = &poly1 * &poly2;
// Dividing polynomials: poly1 / poly2 
let (quotient, rem) = UnivarPolynomial::long_division(&poly1, &poly2);
```

