# amcl_rust_wrapper

- Wraps some parts of [AMCL](https://github.com/miracl/amcl) to provide a nice abstraction to work with finite field elements and group elements when working with elliptic curves.   
- Overloads +, -, *, +=, -= to use with field as well as group elements.  The overloaded operators correspond to constant time methods. But for scalar multiplication, variable time algorithms are present but can be used by calling methods only. 
- Provides abstraction for creating vectors of field elements or group elements and then scale, add, subtract, take inner product or Hadamard product.   
- Additionally, implements some extra algorithms like variable time scalar multiplication using wNAF, constant time and variable time multi-scalar multiplication, batch (simultaneous) inversion and Barrett reduction.

## Building
The wrapper has to be built by enabling any one of the mentioned curve as a feature.   
To build for BLS-381 curve, use

```
cargo build --features bls381
```

Similarly, to build for secp256k1, use
```
cargo build --features secp256k1
```

To run tests for secp256k1, use
```
cargo test --features secp256k1
```

To use it as dependency crate, add the name of the curve as a feature. Something like this
```
[dependencies.amcl_wrapper]
git = "https://github.com/lovesh/amcl_rust_wrapper"
branch = "master"
features = ["bls381"]
```

Note that only one curve can be used at a time so the code only works with one feature.
 
## Examples
1. Create some random field elements or group elements and do some basic additions/subtraction/multiplication 
```
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

```
// Compute inverse of element
let n = a.invert();

// Faster but constant time computation of inverse of several elements. `inverses` will be a vector of inverses of elements and `pr` will be the product of all inverses.
let (inverses, pr) = FieldElement::batch_invert(vec![a, b, c, d].as_slice());  
```

```
// Compute a width 5 NAF.
let a_wnaf = a.to_wnaf(5);
```

```
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

```
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
```

Mutating versions of the above operations like addition/subtraction/negation/inversion are present but have to be called as methods like `b.negate()`

2. Scalar multiplication
```
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
```
// creates a vector of size 10 with all elements as 0
let mut a = FieldElementVector::new(10);
// Add 2 more elements to the above vector
a.push(FieldElement::random());
a.push(FieldElement::random());

a[0];       // 0th element of above vector
a[1];       // 1st element of above vector

a.len();    // length of vector
a.sum();    // sum of elements of vector 
```

```
// Return a Vandermonde vector of a given field element, i.e. given element `k` and size `n`, return vector as `vec![1, k, k^2, k^3, ... k^n-1]`
let k = FieldElement::random();  
let van_vec = FieldElementVector::new_vandermonde_vector(&k, 5);
```

```
// creates a vector of size 10 with randomly generated field elements
let rands: Vec<_> = (0..10).map(|_| FieldElement::random()).collect();

// an alternative way of creating vector of size 10 of random field elements
let rands_1 = FieldElementVector::random(10);
```

```
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

```
// Compute inner product of 2 vectors. Requires vectors to be of equal length.
let ip = rands.inner_product(&rands_1);

// Compute Hadamard product of 2 vectors. Requires vectors to be of equal length.
let hp = rands.hadamard_product(&rands_1);
```

5. Create vectors of group elements and do some operations
```
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

```
// Compute inner product of a vector of group elements with a vector of field elements.    
// eg. given a vector of group elements and field elements, G and F respectively, compute G[0]*F[0] + G[1]*F[1] + G[2]*F[2] + .. G[n-1]*F[n-1]   
// requires vectors to be of same length
let g = G1Vector::random(10);
let f = FieldElementVector::random(10);

// Uses constant time multi-scalar multiplication `multi_scalar_mul_const_time` underneath. 
let ip = g.inner_product_const_time(&f);

// Uses variable time multi-scalar multiplication `multi_scalar_mul_var_time` underneath.
let ip1 = g.inner_product_var_time(&f);

// If lookup tables are already constructed, `multi_scalar_mul_const_time_with_precomputation_done` and `multi_scalar_mul_var_time_with_precomputation_done` can be used for constant and variable time multi-scalar multiplication
```

5. Pairing support. Ate pairing is supported with target group `GT`  
```
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

// Raise target group element to field element (GT^f)
let r = FieldElement::random();
let p = gt.pow(&r);
```