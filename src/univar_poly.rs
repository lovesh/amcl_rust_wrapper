use std::ops::{Add, Index, IndexMut, Mul, Sub};

use crate::field_elem::{FieldElement, FieldElementVector};
use crate::rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use std::cmp::max;

/// Univariate polynomial represented with coefficients in a vector. The ith element of the vector is the coefficient of the ith degree term.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct UnivarPolynomial(pub FieldElementVector);

impl UnivarPolynomial {
    /// Return a zero polynomial of degree `degree`
    pub fn new(degree: usize) -> Self {
        let coeffs = FieldElementVector::new(degree + 1);
        UnivarPolynomial(coeffs)
    }

    /// Return a constant polynomial
    pub fn new_constant(constant: FieldElement) -> Self {
        let mut coeffs = FieldElementVector::new(1);
        coeffs[0] = constant;
        UnivarPolynomial(coeffs)
    }

    /// Return a randomly chosen polynomial (each coefficient is randomly chosen) of degree `degree`.
    pub fn random(degree: usize) -> Self {
        Self(FieldElementVector::random(degree + 1)) // +1 for constant term
    }

    /// Create a polynomial with given roots in `roots`
    /// i.e. (x-roots[0])*(x-roots[1])*(x-roots[2])...(x-roots[last]) given `roots`
    pub fn new_with_roots(roots: &[FieldElement]) -> Self {
        // vector of [(x-roots[0]), (x-roots[1]), (x-roots[2]), ...]
        let x_i = roots
            .par_iter()
            .map(|i| {
                let mut v = FieldElementVector::with_capacity(2);
                v.push(-i);
                v.push(FieldElement::one());
                UnivarPolynomial(v)
            })
            .collect::<Vec<UnivarPolynomial>>();

        // Polynomial (x-roots[0])*(x-roots[1])*(x-roots[2])...(x-roots[last])
        x_i.par_iter().cloned().reduce(
            || Self::new_constant(FieldElement::one()),
            |a, b| UnivarPolynomial::multiply(&a, &b),
        )
    }

    pub fn coefficients(&self) -> &FieldElementVector {
        &self.0
    }

    pub fn degree(&self) -> usize {
        // TODO: This makes fetching the coefficient ambiguous as a 0 degree polynomial might
        // have a coefficient for the 0th degree or it might not. Should probably adapt Index and IndexMut trait.
        let l = self.0.len();
        if l == 0 {
            l
        } else {
            l - 1
        }
    }

    /// Polynomial is zero if all coefficients are 0
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|coeff| coeff.is_zero())
    }

    // Evaluate polynomial at given `x`
    pub fn eval(&self, x: &FieldElement) -> FieldElement {
        if x.is_zero() {
            self[0].clone()
        } else {
            // Use Horner's method https://en.wikipedia.org/wiki/Horner%27s_method
            // p(x) = a_0 + a_1*x + a_2*x^2 + a_3*x^3 + a_4*x^4 + ...
            // p(x) = a_0 + x*(a_1 + x*(a_2 + x*(a_3 + x*(a_4 + ... x*(a_{n-1} + x*a_n))))..
            // Reading coefficients from higher to lower degrees.
            let mut res = self.0[self.0.len() - 1].clone(); // a_n
            for i in (0..=self.0.len() - 2).rev() {
                // in each iteration, multiply `res` with `x` and add the coefficient for ith degree, a_i
                res = &self.0[i] + &(&res * x);
            }
            res
        }
    }

    /// Divides 2 polynomials i.e. `dividend` / `divisor` using long division.
    /// Returns (quotient, remainder)
    pub fn long_division(dividend: &Self, divisor: &Self) -> (Self, Self) {
        assert!(!divisor.is_zero());
        assert!(!divisor[divisor.degree()].is_zero());

        let mut remainder: UnivarPolynomial = dividend.clone();
        let mut quotient = vec![];
        // Inverse of coefficient of highest degree of the divisor polynomial. This will be multiplied
        // with the coefficient of highest degree of the remainder.
        let highest_degree_coeff_inv = divisor[divisor.degree()].inverse();
        let rem_degree = dividend.degree();
        let div_degree = divisor.degree();
        for i in (div_degree..=rem_degree).rev() {
            if remainder[i].is_zero() {
                quotient.push(FieldElement::zero());
                continue;
            }

            let q = &highest_degree_coeff_inv * &remainder[i];
            for j in 0..div_degree {
                remainder[i - div_degree + j] -= &(&divisor[j] * &q);
            }
            quotient.push(q);
        }
        // The coefficients of the quotient polynomial were computed from highest to lowest degree.
        quotient.reverse();
        // Remainder's degree will be less than divisor's degree.
        for _ in div_degree..=rem_degree {
            remainder.0.pop();
        }
        (
            UnivarPolynomial(FieldElementVector::from(quotient)),
            remainder,
        )
    }

    /// Return product of 2 polynomials. `left` * `right`
    pub fn multiply(left: &Self, right: &Self) -> Self {
        let mut product = Self::new(left.degree() + right.degree());
        for i in 0..=left.degree() {
            for j in 0..=right.degree() {
                product[i + j] += &left[i] * &right[j];
            }
        }
        product
    }

    /// Return sum of 2 polynomials. `left` + `right`
    pub fn sum(left: &Self, right: &Self) -> Self {
        // The resulting sum polynomial is initialized with the input polynomial of larger degree
        let (mut sum_poly, smaller_poly, smaller_poly_degree) = if left.degree() > right.degree() {
            (left.clone(), right, right.degree())
        } else {
            (right.clone(), left, left.degree())
        };

        // The following unobvious code is to use rayon for parallelization. A simpler (non-parallel)
        // version would be  `for i in 0..=smaller_poly_degree { sum_poly[i] += &smaller_poly[i]; }`

        // Add small degree ([0, smaller_poly_degree]) terms in parallel
        let small_degree_terms = (0..=smaller_poly_degree)
            .into_par_iter()
            .map(|i| &sum_poly[i] + &smaller_poly[i])
            .collect::<Vec<FieldElement>>();
        // Replace small degree ([0, smaller_poly_degree]) terms in the sum_poly
        sum_poly.replace_small_degree_terms(smaller_poly_degree, small_degree_terms.into_iter());
        sum_poly
    }

    /// Return difference of 2 polynomials. `left` - `right`
    pub fn difference(left: &Self, right: &Self) -> Self {
        let left_degree = left.degree();
        let right_degree = right.degree();
        let diff_poly_degree = max(left_degree, right_degree);
        let mut diff = Self::new(diff_poly_degree);
        for i in 0..=diff_poly_degree {
            if i <= left_degree {
                diff[i] = left[i].clone();
            }
            if i <= right_degree {
                diff[i] -= &right[i];
            }
        }
        diff
    }

    pub fn multiply_by_constant(&self, constant: &FieldElement) -> UnivarPolynomial {
        let mut new_poly = self.clone();
        for i in 0..=self.degree() {
            new_poly[i] = constant * &self[i];
        }
        new_poly
    }

    pub fn multiply_by_monic_monomial(&self, monomial_degree: u64) -> UnivarPolynomial {
        let mut new_poly = self.clone();
        let new_poly_beginning = FieldElementVector::new(monomial_degree as usize);
        new_poly.0.splice(0..0, new_poly_beginning);
        new_poly
    }

    /// Replace terms of `self` from degree 0 to `till_degree` with coefficients in `replace_with`.
    /// Assumes `replace_with` will yield at least `till_degree` + 1 coefficients
    fn replace_small_degree_terms<I: IntoIterator<Item = FieldElement>>(
        &mut self,
        till_degree: usize,
        replace_with: I,
    ) {
        self.0.splice(0..=till_degree, replace_with)
    }
}

impl Index<usize> for UnivarPolynomial {
    type Output = FieldElement;

    fn index(&self, idx: usize) -> &FieldElement {
        &self.0[idx]
    }
}

impl IndexMut<usize> for UnivarPolynomial {
    fn index_mut(&mut self, idx: usize) -> &mut FieldElement {
        &mut self.0[idx]
    }
}

impl Eq for UnivarPolynomial {}

impl<'a> Add<&'a UnivarPolynomial> for &UnivarPolynomial {
    type Output = UnivarPolynomial;
    fn add(self, other: &'a UnivarPolynomial) -> UnivarPolynomial {
        UnivarPolynomial::sum(self, other)
    }
}

impl<'a> Sub<&'a UnivarPolynomial> for &UnivarPolynomial {
    type Output = UnivarPolynomial;

    fn sub(self, other: &'a UnivarPolynomial) -> UnivarPolynomial {
        UnivarPolynomial::difference(self, other)
    }
}

impl<'a> Mul<&'a UnivarPolynomial> for &UnivarPolynomial {
    type Output = UnivarPolynomial;

    fn mul(self, other: &'a UnivarPolynomial) -> UnivarPolynomial {
        UnivarPolynomial::multiply(self, other)
    }
}

/// Creates a new univariate polynomial from given coefficients from lower to higher degree terms
#[macro_export]
macro_rules! univar_polynomial {
    ( $( $elem:expr ),* ) => {
        {
            let mut coeffs = vec![];
            $(
                coeffs.push($elem);
            )*
            UnivarPolynomial(coeffs.into())
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use std::time::{Duration, Instant};

    #[test]
    fn test_poly() {
        let degree = 10;
        let poly1 = UnivarPolynomial(FieldElementVector::random(degree + 1));
        assert!(!poly1.is_zero());

        let poly2 = UnivarPolynomial(FieldElementVector::new(degree + 1));
        assert!(poly2.is_zero());

        let poly3 = UnivarPolynomial::new(degree);
        assert!(poly3.is_zero());

        let poly4 = UnivarPolynomial::new_constant(FieldElement::from(100u64));
        assert!(!poly4.is_zero());
        assert_eq!(poly4.degree(), 0);
        assert_eq!(poly4[0], FieldElement::from(100u64));
    }

    #[test]
    fn test_create_poly_from_macro() {
        let poly = univar_polynomial!(
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::from(87u64),
            FieldElement::minus_one(),
            FieldElement::from(300u64)
        );
        assert_eq!(poly.degree(), 4);
        assert_eq!(poly[0], FieldElement::one());
        assert_eq!(poly[1], FieldElement::zero());
        assert_eq!(poly[2], FieldElement::from(87u64));
        assert_eq!(poly[3], FieldElement::minus_one());
        assert_eq!(poly[4], FieldElement::from(300u64));
    }

    #[test]
    fn test_poly_long_div() {
        // x^2 - 1 / x + 1 = x - 1
        // dividend = -1 + x^2
        let c1 = vec![
            FieldElement::minus_one(),
            FieldElement::zero(),
            FieldElement::one(),
        ];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = 1 + x
        let c2 = vec![FieldElement::one(), FieldElement::one()];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, rem) = UnivarPolynomial::long_division(&dividend, &divisor);
        println!("Quotient={:?}", &quotient);
        // quotient = -1 + x
        assert_eq!(quotient.degree(), 1);
        assert_eq!(quotient[0], FieldElement::minus_one());
        assert_eq!(quotient[1], FieldElement::one());

        assert_eq!(rem.degree(), 0);

        let quotient = UnivarPolynomial::long_division(&dividend, &quotient).0;
        println!("Quotient={:?}", &quotient);
        // quotient = 1 + x
        assert_eq!(quotient.degree(), 1);
        assert_eq!(quotient[0], FieldElement::one());
        assert_eq!(quotient[1], FieldElement::one());

        // 2x^2 + 3x + 1 / x + 1 = 2x + 1
        // dividend = 1 + 3x + 2x^2
        let c1 = vec![
            FieldElement::one(),
            FieldElement::from(3u64),
            FieldElement::from(2u64),
        ];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = 1 + x
        let c2 = vec![FieldElement::one(), FieldElement::one()];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, rem) = UnivarPolynomial::long_division(&dividend, &divisor);
        println!("Quotient={:?}", &quotient);
        // quotient = 1 + 2x
        assert_eq!(quotient.degree(), 1);
        assert_eq!(quotient[0], FieldElement::one());
        assert_eq!(quotient[1], FieldElement::from(2u64));

        assert_eq!(rem.degree(), 0);

        // 4x - 4 / x - 1 = 4
        // dividend = -4 + 4x
        let c1 = vec![-FieldElement::from(4u64), FieldElement::from(4u64)];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = -1 + x
        let c2 = vec![FieldElement::minus_one(), FieldElement::one()];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, rem) = UnivarPolynomial::long_division(&dividend, &divisor);
        println!("Quotient={:?}", &quotient);

        // quotient = 4
        assert_eq!(quotient.degree(), 0);
        assert_eq!(quotient[0], FieldElement::from(4u64));

        assert_eq!(rem.degree(), 0);

        // x^5 + x^3 + 4x^2 + 4 / x^2 + 1 = x^3 + 4
        // dividend = 4 + 4x^2 + x^3 + x^5
        let c1 = vec![
            FieldElement::from(4u64),
            FieldElement::zero(),
            FieldElement::from(4u64),
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::one(),
        ];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = 1 + x^2
        let c2 = vec![
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::one(),
        ];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, rem) = UnivarPolynomial::long_division(&dividend, &divisor);
        println!("Quotient={:?}", &quotient);

        // quotient = 4 + x^3
        assert_eq!(quotient.degree(), 3);
        assert_eq!(quotient[0], FieldElement::from(4u64));
        assert_eq!(quotient[1], FieldElement::zero());
        assert_eq!(quotient[2], FieldElement::zero());
        assert_eq!(quotient[3], FieldElement::one());

        assert_eq!(rem.degree(), 1);

        // 2x^4 - 40x^3 + 3x^2 - 56x - 80 / x - 20 = 2x^3 + 3x + 4
        // dividend = -80 - 56x + 3x^2 - 40x^3 + 2x^4
        let c1 = vec![
            -FieldElement::from(80u64),
            -FieldElement::from(56u64),
            FieldElement::from(3u64),
            -FieldElement::from(40u64),
            FieldElement::from(2u64),
        ];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = -20 + x
        let c2 = vec![-FieldElement::from(20), FieldElement::one()];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, rem) = UnivarPolynomial::long_division(&dividend, &divisor);
        println!("Quotient={:?}", &quotient);

        // quotient = 4 + 3x + 2x^3
        assert_eq!(quotient.degree(), 3);
        assert_eq!(quotient[0], FieldElement::from(4u64));
        assert_eq!(quotient[1], FieldElement::from(3u64));
        assert_eq!(quotient[2], FieldElement::zero());
        assert_eq!(quotient[3], FieldElement::from(2u64));

        assert_eq!(rem.degree(), 0);
    }

    #[test]
    fn test_poly_multiply() {
        // (x + 1) * (x - 1) = x^2 - 1
        // x + 1
        let left = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::one(),
            FieldElement::one(),
        ]));
        // -1 + x
        let right = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::minus_one(),
            FieldElement::one(),
        ]));
        let product = UnivarPolynomial::multiply(&left, &right);
        // product = -1 + x^2
        assert_eq!(product.degree(), 2);
        assert_eq!(product[0], FieldElement::minus_one());
        assert_eq!(product[1], FieldElement::zero());
        assert_eq!(product[2], FieldElement::one());

        // Test overloaded operator
        assert_eq!(product, &left * &right);

        // (x + 1) * (2x + 1) = 2x^2 + 3x + 1
        // 1 + x
        let left = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::one(),
            FieldElement::one(),
        ]));
        // 1 + 2x
        let right = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::one(),
            FieldElement::from(2u64),
        ]));
        let product = UnivarPolynomial::multiply(&left, &right);
        // product = 2x^2 + 3x + 1
        assert_eq!(product.degree(), 2);
        assert_eq!(product[0], FieldElement::one());
        assert_eq!(product[1], FieldElement::from(3u64));
        assert_eq!(product[2], FieldElement::from(2u64));

        // Test overloaded operator
        assert_eq!(product, &left * &right);

        // (x^2 + 1) * (x^3 + 4) = x^5 + x^3 + 4x^2 + 4
        // 1 + x^2
        let left = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::one(),
        ]));
        // 4 + x^3
        let right = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::from(4u64),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::one(),
        ]));
        let product = UnivarPolynomial::multiply(&left, &right);
        // 4 + 4x^2 + x^3 + x^5
        assert_eq!(product.degree(), 5);
        assert_eq!(product[0], FieldElement::from(4u64));
        assert_eq!(product[1], FieldElement::zero());
        assert_eq!(product[2], FieldElement::from(4u64));
        assert_eq!(product[3], FieldElement::one());
        assert_eq!(product[4], FieldElement::zero());
        assert_eq!(product[5], FieldElement::one());

        // Test overloaded operator
        assert_eq!(product, &left * &right);
    }

    #[test]
    fn test_poly_rem() {
        // x^2 - 5 / x + 1 => q = x - 1, r = -4
        // dividend = -5 + x^2
        let c1 = vec![
            -FieldElement::from(5u64),
            FieldElement::zero(),
            FieldElement::one(),
        ];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = 1 + x
        let c2 = vec![FieldElement::one(), FieldElement::one()];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, remainder) = UnivarPolynomial::long_division(&dividend, &divisor);
        // quotient = -1 + x
        assert_eq!(quotient.degree(), 1);
        assert_eq!(quotient[0], FieldElement::minus_one());
        assert_eq!(quotient[1], FieldElement::one());

        // remainder = -4
        assert_eq!(remainder.degree(), 0);
        assert_eq!(remainder[0], -FieldElement::from(4u64));

        // x^5 + 2x^3 + 4x^2 + 4 / x^2 + 1 = q = x^3 + x + 4, r = -x
        // dividend = 4 + 4x^2 + 2x^3 + x^5
        let c1 = vec![
            FieldElement::from(4u64),
            FieldElement::zero(),
            FieldElement::from(4u64),
            FieldElement::from(2u64),
            FieldElement::zero(),
            FieldElement::one(),
        ];
        let dividend = UnivarPolynomial(FieldElementVector::from(c1));
        // divisor = 1 + x^2
        let c2 = vec![
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::one(),
        ];
        let divisor = UnivarPolynomial(FieldElementVector::from(c2));
        let (quotient, remainder) = UnivarPolynomial::long_division(&dividend, &divisor);

        // quotient = 4 + x^3
        assert_eq!(quotient.degree(), 3);
        assert_eq!(quotient[0], FieldElement::from(4u64));
        assert_eq!(quotient[1], FieldElement::one());
        assert_eq!(quotient[2], FieldElement::zero());
        assert_eq!(quotient[3], FieldElement::one());

        assert_eq!(remainder.degree(), 1);
        assert_eq!(remainder[0], FieldElement::zero());
        assert_eq!(remainder[1], FieldElement::minus_one());
    }

    #[test]
    fn test_random_poly_sum_difference() {
        // Test sum and difference of randomly generated polynomials.
        let num_test_cases = 100;
        let mut rng = rand::thread_rng();
        let mut start = Instant::now();
        for _ in 0..num_test_cases {
            let left = UnivarPolynomial::random(rng.gen_range(1, 100));
            let right = UnivarPolynomial::random(rng.gen_range(1, 100));
            let sum = UnivarPolynomial::sum(&left, &right);

            // sum is commutative
            assert_eq!(sum, UnivarPolynomial::sum(&right, &left));

            // Test overloaded operator
            assert_eq!(sum, &left + &right);

            // sum - left == right
            let mut diff_1 = UnivarPolynomial::difference(&sum, &right);

            // Test overloaded operator
            assert_eq!(diff_1, &sum - &right);

            // Since degree of difference is same as degree of `sum` but the higher degree coeffs
            // of difference will be 0. Remove those 0s (after checking that they really are 0) and
            // then do equality comparison with `left`
            while diff_1.degree() > left.degree() {
                let c = diff_1.0.pop().unwrap();
                assert!(c.is_zero());
            }
            assert_eq!(diff_1, left);

            // sum - right == left
            let mut diff_2 = UnivarPolynomial::difference(&sum, &left);

            // Test overloaded operator
            assert_eq!(diff_2, &sum - &left);

            // Since degree of difference is same as degree of `sum` but the higher degree coeffs
            // of difference will be 0. Remove those 0s (after checking that they really are 0) and
            // then do equality comparison with `right`
            while diff_2.degree() > right.degree() {
                let c = diff_2.0.pop().unwrap();
                assert!(c.is_zero());
            }
            assert_eq!(diff_2, right);
        }
        println!(
            "Sum diff time for {} elems = {:?}",
            num_test_cases,
            start.elapsed()
        );
    }

    #[test]
    fn test_random_poly_long_div() {
        // Multiply 2 random polynomials and then use the result to check long division
        let num_test_cases = 100;
        let mut rng = rand::thread_rng();
        for _ in 0..num_test_cases {
            let left = UnivarPolynomial::random(rng.gen_range(1, 100));
            let right = UnivarPolynomial::random(rng.gen_range(1, 100));
            let product = UnivarPolynomial::multiply(&left, &right);

            // product / left == right
            let quotient_1 = UnivarPolynomial::long_division(&product, &left).0;
            assert_eq!(quotient_1, right);

            // product / right == left
            let quotient_2 = UnivarPolynomial::long_division(&product, &right).0;
            assert_eq!(quotient_2, left);

            // Test overloaded operator
            assert_eq!(product, &left * &right);
        }
    }

    #[test]
    fn test_random_poly_long_div_remainder() {
        // Divide 2 random polynomials and check that the quotient and remainder are correct using
        // the relation dividend = divisor * quotient + remainder
        let num_test_cases = 100;
        let mut rng = rand::thread_rng();
        for _ in 0..num_test_cases {
            let d_1: usize = rng.gen_range(1, 100);
            let d_2: usize = rng.gen_range(1, 100);
            let (dividend, divisor) = if d_1 > d_2 {
                (UnivarPolynomial::random(d_1), UnivarPolynomial::random(d_2))
            } else {
                (UnivarPolynomial::random(d_2), UnivarPolynomial::random(d_1))
            };
            // dividend / divisor => quotient and remainder
            let (quotient, remainder) = UnivarPolynomial::long_division(&dividend, &divisor);

            // dividend = divisor * quotient + remainder

            // div_quo = divisor * quotient
            let div_quo = UnivarPolynomial::multiply(&divisor, &quotient);
            // expected_dividend = div_quo + remainder
            let expected_dividend = UnivarPolynomial::sum(&div_quo, &remainder);
            assert_eq!(expected_dividend, dividend);
        }
    }

    #[test]
    fn test_poly_from_given_roots() {
        // Check resulting polynomial is of correct degree and polynomial becomes 0 at each root
        let num_test_cases = 100;
        let mut rng = rand::thread_rng();
        let mut start = Instant::now();
        for _ in 0..num_test_cases {
            let num_roots = rng.gen_range(2, 30);
            let roots = FieldElementVector::random(num_roots);
            let poly = UnivarPolynomial::new_with_roots(roots.as_slice());
            assert_eq!(poly.degree(), num_roots);
            for r in roots {
                assert_eq!(poly.eval(&r), FieldElement::zero())
            }
        }
        println!("Time for {} elems = {:?}", num_test_cases, start.elapsed());
    }

    #[test]
    fn test_multiply_with_constant() {
        // 9 + 2x + 75x^2 + 128x^3
        let orig = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::from(9u64),
            FieldElement::from(2u64),
            FieldElement::from(75u64),
            FieldElement::from(128u64),
        ]));
        let c = FieldElement::from(3u64);
        let new = orig.multiply_by_constant(&c);
        assert_eq!(new.degree(), 3);
        assert_eq!(new[0], FieldElement::from(27));
        assert_eq!(new[1], FieldElement::from(6));
        assert_eq!(new[2], FieldElement::from(225));
        assert_eq!(new[3], FieldElement::from(384));

        // 1 + 4x^2 + 5x^3 + 18x^6
        let orig = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::from(4u64),
            FieldElement::from(5u64),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::from(18u64),
        ]));
        let c = FieldElement::from(10u64);
        let new = orig.multiply_by_constant(&c);
        assert_eq!(new.degree(), 6);
        assert_eq!(new[0], FieldElement::from(10));
        assert_eq!(new[1], FieldElement::zero());
        assert_eq!(new[2], FieldElement::from(40));
        assert_eq!(new[3], FieldElement::from(50));
        assert_eq!(new[4], FieldElement::zero());
        assert_eq!(new[5], FieldElement::zero());
        assert_eq!(new[6], FieldElement::from(180));

        // take a random polynomial, multiply it with a constant, then multiply it with inverse of
        // the same constant. result should be same as original
        let random_poly = UnivarPolynomial::random(10);
        let c = FieldElement::random();
        let c_inv = c.inverse();
        assert_eq!(
            random_poly,
            random_poly
                .multiply_by_constant(&c)
                .multiply_by_constant(&c_inv)
        );
    }

    #[test]
    fn test_multiply_with_monic_monomial() {
        // 9 + 2x + 75x^2 + 128x^3
        let orig = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::from(9u64),
            FieldElement::from(2u64),
            FieldElement::from(75u64),
            FieldElement::from(128u64),
        ]));

        let monomial_degree = 0;
        let new = orig.multiply_by_monic_monomial(monomial_degree);
        assert_eq!(new, orig);

        let monomial_degree = 1;
        let new = orig.multiply_by_monic_monomial(monomial_degree);
        assert_eq!(new.degree(), 4);
        assert_eq!(new[0], FieldElement::zero());
        assert_eq!(new[1], FieldElement::from(9u64));
        assert_eq!(new[2], FieldElement::from(2u64));
        assert_eq!(new[3], FieldElement::from(75u64));
        assert_eq!(new[4], FieldElement::from(128u64));

        let monomial_degree = 2;
        let new = orig.multiply_by_monic_monomial(monomial_degree);
        assert_eq!(new.degree(), 5);
        assert_eq!(new[0], FieldElement::zero());
        assert_eq!(new[1], FieldElement::zero());
        assert_eq!(new[2], FieldElement::from(9u64));
        assert_eq!(new[3], FieldElement::from(2u64));
        assert_eq!(new[4], FieldElement::from(75u64));
        assert_eq!(new[5], FieldElement::from(128u64));

        // 1 + 4x^2 + 5x^3 + 18x^6
        let orig = UnivarPolynomial(FieldElementVector::from(vec![
            FieldElement::one(),
            FieldElement::zero(),
            FieldElement::from(4u64),
            FieldElement::from(5u64),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::from(18u64),
        ]));

        let monomial_degree = 0;
        let new = orig.multiply_by_monic_monomial(monomial_degree);
        assert_eq!(new, orig);

        let monomial_degree = 1;
        let new = orig.multiply_by_monic_monomial(monomial_degree);
        assert_eq!(new.degree(), 7);
        assert_eq!(new[0], FieldElement::zero());
        assert_eq!(new[1], FieldElement::one());
        assert_eq!(new[2], FieldElement::zero());
        assert_eq!(new[3], FieldElement::from(4));
        assert_eq!(new[4], FieldElement::from(5));
        assert_eq!(new[5], FieldElement::zero());
        assert_eq!(new[6], FieldElement::zero());
        assert_eq!(new[7], FieldElement::from(18));

        let monomial_degree = 2;
        let new = orig.multiply_by_monic_monomial(monomial_degree);
        assert_eq!(new.degree(), 8);
        assert_eq!(new[0], FieldElement::zero());
        assert_eq!(new[1], FieldElement::zero());
        assert_eq!(new[2], FieldElement::one());
        assert_eq!(new[3], FieldElement::zero());
        assert_eq!(new[4], FieldElement::from(4));
        assert_eq!(new[5], FieldElement::from(5));
        assert_eq!(new[6], FieldElement::zero());
        assert_eq!(new[7], FieldElement::zero());
        assert_eq!(new[8], FieldElement::from(18));
    }
}
