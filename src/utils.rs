extern crate rand;
extern crate sha3;

use rand::prelude::*;

use crate::constants::{CURVE_ORDER, MODBYTES};
use crate::types::{BigNum, DoubleBigNum};

use sha3::Shake256;
use sha3::digest::{Input, ExtendableOutput, XofReader};


/// Generate a random number twice the number of bytes as the curve reduced by mod CURVE_ORDER
/// Here 2X bytes are randomly generated mod CURVE_ORDER. The result is CURVE_ORDER bytes.
pub fn random_mod_order<R: RngCore>(r: Option<&mut R>) -> BigNum {
    let mut seed1 = vec![0u8; MODBYTES];
    let mut seed2 = vec![0u8; MODBYTES];

    match r {
        Some(rr) => {
            rr.fill_bytes(&mut seed1.as_mut_slice());
            rr.fill_bytes(&mut seed2.as_mut_slice());
        }
        None => {
            thread_rng().fill_bytes(&mut seed1.as_mut_slice());
            thread_rng().fill_bytes(&mut seed2.as_mut_slice());
        }
    }

    compute_big(seed1.as_slice(), seed2.as_slice())
}

/// Hash the data to a field order element
/// Based on standardization efforts for BLS12-381, it is recommended to use XOF
/// instead of SHA2 due to its properties of better hiding low entropy inputs like HKDF or SHAKE.
///
/// Here, 2X bytes is used. This allows a bigger range than SHA2-256 and eliminates some biases
/// for outputs that are close to the CURVE_ORDER like SHA2-384.
///
/// Salt should be of sufficient length to add the security of the keying material (See https://tools.ietf.org/html/rfc5869)
/// The output of 2X bytes should require the salt to not be less than 32 bytes
///
/// The result is CURVE_ORDER bytes.
pub fn hash_mod_order(data: &[u8]) -> BigNum {
    let mut hasher = Shake256::default();
    hasher.input(data);
    let mut output = vec![0u8; MODBYTES * 2];
    hasher.xof_result().read(&mut output);

    compute_big(
        &output.as_slice()[0..MODBYTES],
        &output.as_slice()[MODBYTES..],
    )
}

fn compute_big(seed1: &[u8], seed2: &[u8]) -> BigNum {
    let num1 = BigNum::frombytes(seed1);
    let num2 = BigNum::frombytes(seed2);
    let num1 = DoubleBigNum::new_scopy(&num1);
    let mut res = DoubleBigNum::new();
    res.ucopy(&num2);
    res.add(&num1);
    res.dmod(&CURVE_ORDER)
}


/// Hash message and return output of size equal to curve modulus. Uses SHAKE to hash the message.
pub fn hash_msg(msg: &[u8]) -> [u8; MODBYTES] {
    let mut hasher = Shake256::default();
    hasher.input(&msg);
    let mut h: [u8; MODBYTES] = [0; MODBYTES];
    hasher.xof_result().read(&mut h);
    h
}

/// Perform Barrett reduction given the params computed from `barrett_reduction_params`. Algorithm 14.42 from Handbook of Applied Cryptography
pub fn barrett_reduction(
    x: &DoubleBigNum,
    modulus: &BigNum,
    k: usize,
    u: &BigNum,
    v: &BigNum,
) -> BigNum {
    // q1 = floor(x / 2^{k-1})
    let mut q1 = x.clone();
    q1.shr(k - 1);
    // Above right shift will convert q from DBIG to BIG
    let q1 = BigNum::new_dcopy(&q1);

    let q2 = BigNum::mul(&q1, &u);

    // q3 = floor(q2 / 2^{k+1})
    let mut q3 = q2.clone();
    q3.shr(k + 1);
    let q3 = BigNum::new_dcopy(&q3);

    // r1 = x % 2^{k+1}
    let mut r1 = x.clone();
    r1.mod2m(k + 1);
    let r1 = BigNum::new_dcopy(&r1);

    // r2 = (q3 * modulus) % 2^{k+1}
    let mut r2 = BigNum::mul(&q3, modulus);
    r2.mod2m(k + 1);
    let r2 = BigNum::new_dcopy(&r2);

    // if r1 > r2, r = r1 - r2 else r = r1 - r2 + v
    // Since negative numbers are not supported, use r2 - r1. This holds since r = r1 - r2 + v = v - (r2 - r1)
    let diff = BigNum::comp(&r1, &r2);
    //println!("diff={}", &diff);
    let mut r = if diff < 0 {
        let m = r2.minus(&r1);
        v.minus(&m)
    } else {
        r1.minus(&r2)
    };
    r.norm();

    // while r >= modulus, r = r - modulus
    while BigNum::comp(&r, modulus) >= 0 {
        r = BigNum::minus(&r, modulus);
        r.norm();
    }
    r
}

// Reducing BigNum for comparison with `rmod`
fn __barrett_reduction__(x: &BigNum, modulus: &BigNum, k: usize, u: &BigNum, v: &BigNum) -> BigNum {
    // q1 = floor(x / 2^{k-1})
    let mut q1 = x.clone();
    q1.shr(k - 1);

    let q2 = BigNum::mul(&q1, &u);

    // q3 = floor(q2 / 2^{k+1})
    let mut q3 = q2.clone();
    q3.shr(k + 1);
    let q3 = BigNum::new_dcopy(&q3);

    // r1 = x % 2^{k+1}
    let mut r1 = x.clone();
    r1.mod2m(k + 1);

    // r2 = (q3 * modulus) % 2^{k+1}
    let mut r2 = BigNum::mul(&q3, modulus);
    r2.mod2m(k + 1);
    let r2 = BigNum::new_dcopy(&r2);

    // if r1 > r2, r = r1 - r2 else r = r1 - r2 + v
    // Since negative numbers are not supported, use r2 - r1. This holds since r = r1 - r2 + v = v - (r2 - r1)
    let diff = BigNum::comp(&r1, &r2);
    //println!("diff={}", &diff);
    let mut r = if diff < 0 {
        let m = r2.minus(&r1);
        v.minus(&m)
    } else {
        r1.minus(&r2)
    };
    r.norm();

    // while r >= modulus, r = r - modulus
    while BigNum::comp(&r, modulus) >= 0 {
        r = BigNum::minus(&r, modulus);
        r.norm();
    }
    r
}

/// For a modulus returns
/// k = number of bits in modulus
/// u = floor(2^2k / modulus)
/// v = 2^(k+1)
pub fn barrett_reduction_params(modulus: &BigNum) -> (usize, BigNum, BigNum) {
    let k = modulus.nbits();

    // u = floor(2^2k/CURVE_ORDER)
    let mut u = DoubleBigNum::new();
    u.w[0] = 1;
    // `u.shl(2*k)` crashes, so perform shl(k) twice
    u.shl(k);
    u.shl(k);
    // div returns floored value
    let u = u.div(&CURVE_ORDER);

    // v = 2^(k+1)
    let mut v = BigNum::new_int(1isize);
    v.shl(k + 1);

    (k, u, v)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::constants;
    use crate::field_elem::FieldElement;
    use crate::group_elem::GroupElement;
    use crate::group_elem_g1::G1;
    use crate::utils::rand::Rng;
    use crate::ECCurve::big::BIG;
    use crate::ECCurve::ecp::ECP;
    use crate::ECCurve::fp::FP;
    use rand_chacha::ChaChaRng;
    use rand_core::SeedableRng;
    use std::time::Instant;

    #[test]
    fn timing_fp_big() {
        // TODO: Compare adding raw BIGs and FieldElement to check the overhead of the abstraction
        let count = 100;
        let elems: Vec<_> = (0..count).map(|_| FieldElement::random()).collect();
        let bigs: Vec<_> = elems.iter().map(|f| f.to_bignum()).collect();
        let fs: Vec<_> = bigs.iter().map(|b| FP::new_big(&b)).collect();
        let mut res_mul = BIG::new_int(1);
        let mut start = Instant::now();
        for b in &bigs {
            res_mul = BigNum::modmul(&res_mul, &b, &CURVE_ORDER);
        }
        println!(
            "Multiplication time for {} BIGs = {:?}",
            count,
            start.elapsed()
        );

        let mut res_mul = FP::new_int(1);
        start = Instant::now();
        for f in &fs {
            res_mul.mul(&f);
        }
        println!(
            "Multiplication time for {} FPs = {:?}",
            count,
            start.elapsed()
        );

        let res_mul = FieldElement::one();
        start = Instant::now();
        for e in &elems {
            res_mul.multiply(&e);
        }
        println!(
            "Multiplication time for {} FieldElements = {:?}",
            count,
            start.elapsed()
        );

        let mut inverses_b: Vec<BigNum> = vec![];
        let mut inverses_f: Vec<FP> = vec![];

        start = Instant::now();
        for b in &bigs {
            let mut i = b.clone();
            i.invmodp(&CURVE_ORDER);
            inverses_b.push(i);
        }
        println!("Inverse time for {} BIGs = {:?}", count, start.elapsed());
        for i in 0..count {
            let r = BigNum::modmul(&inverses_b[i], &bigs[i], &CURVE_ORDER);
            assert_eq!(BigNum::comp(&r, &BigNum::new_int(1)), 0);
        }

        start = Instant::now();
        for f in &fs {
            let mut i = f.clone();
            i.inverse();
            inverses_f.push(i);
        }
        println!("Inverse time for {} FPs = {:?}", count, start.elapsed());
        for i in 0..count {
            let mut c = inverses_f[i].clone();
            c.mul(&fs[i]);
            assert!(c.equals(&FP::new_int(1)));
        }

        // Fixme: add in FP crashes while adding 100 elems
        let c = 50;
        start = Instant::now();
        let mut r = bigs[0];
        for i in 0..c {
            r.add(&bigs[i]);
            r.rmod(&CURVE_ORDER);
        }
        println!("Addition time for {} BIGs = {:?}", c, start.elapsed());

        let mut r1 = fs[0];
        start = Instant::now();
        for i in 0..c {
            r1.add(&fs[i]);
        }
        println!("Addition time for {} FPs = {:?}", c, start.elapsed());
    }

    #[test]
    fn timing_ecp() {
        let count = 100;
        let mut a = vec![];
        let mut b = vec![];
        let mut g = Vec::<ECP>::new();
        let mut h = Vec::<ECP>::new();

        let mut r1 = vec![];
        let mut r2 = vec![];

        for _ in 0..count {
            a.push(FieldElement::random().to_bignum());
            b.push(FieldElement::random().to_bignum());
            let mut x: G1 = GroupElement::random();
            g.push(x.to_ecp());
            x = GroupElement::random();
            h.push(x.to_ecp());
        }

        let mut start = Instant::now();
        for i in 0..count {
            r1.push(g[i].mul2(&a[i], &h[i], &b[i]));
        }
        println!("mul2 time for {} = {:?}", count, start.elapsed());

        start = Instant::now();
        for i in 0..count {
            let mut _1 = g[i].mul(&a[i]);
            _1.add(&h[i].mul(&b[i]));
            r2.push(_1);
        }
        println!("mul+add time for {} = {:?}", count, start.elapsed());

        for i in 0..count {
            assert!(r1[i].equals(&mut r2[i]))
        }
    }

    #[test]
    fn timing_barrett_reduction() {
        //let (k, u, v) = barrett_reduction_params(&CURVE_ORDER);
        let (k, u, v) = (
            *constants::BARRETT_REDC_K,
            *constants::BARRETT_REDC_U,
            *constants::BARRETT_REDC_V,
        );
        let mut xs = vec![];
        let mut reduced1 = vec![];
        let mut reduced2 = vec![];
        let mut rng = thread_rng();
        let count = 1000;
        for _ in 0..count {
            let a: u32 = rng.gen();
            let s = BigNum::new_int(a as isize);
            let _x = CURVE_ORDER.minus(&s);
            xs.push(BigNum::mul(&_x, &_x));
        }

        let mut start = Instant::now();
        for x in &xs {
            let r = barrett_reduction(&x, &CURVE_ORDER, k, &u, &v);
            reduced1.push(r);
        }
        println!("Barrett time = {:?}", start.elapsed());

        start = Instant::now();
        for x in &xs {
            let mut y = x.clone();
            let z = y.dmod(&CURVE_ORDER);
            reduced2.push(z);
        }
        println!("Normal time = {:?}", start.elapsed());

        for i in 0..count {
            assert_eq!(BigNum::comp(&reduced1[i], &reduced2[i]), 0);
        }
    }

    #[test]
    fn timing_rmod_with_barrett_reduction() {
        let (k, u, v) = (
            *constants::BARRETT_REDC_K,
            *constants::BARRETT_REDC_U,
            *constants::BARRETT_REDC_V,
        );
        let count = 100;
        let elems: Vec<_> = (0..count).map(|_| FieldElement::random()).collect();
        let bigs: Vec<_> = elems.iter().map(|f| f.to_bignum()).collect();

        let mut sum = bigs[0].clone();
        let mut start = Instant::now();
        for i in 0..count {
            sum = BigNum::plus(&sum, &bigs[i]);
            sum.rmod(&CURVE_ORDER)
        }
        println!("rmod time = {:?}", start.elapsed());

        let mut sum_b = bigs[0].clone();
        start = Instant::now();
        for i in 0..count {
            sum_b = BigNum::plus(&sum_b, &bigs[i]);
            sum_b = __barrett_reduction__(&sum_b, &CURVE_ORDER, k, &u, &v)
        }
        println!("Barrett time = {:?}", start.elapsed());

        assert_eq!(BigNum::comp(&sum, &sum_b), 0)
    }

    #[test]
    fn random_tests() {
        let seed = b"11111111111111111111111111111111";
        let mut rng = ChaChaRng::from_seed(*seed);

        let r = random_mod_order(Some(&mut rng));
        assert_eq!(r, BigNum::from_hex("000000000000000000000000000000005C9F002063FDC7EDD33E8787C8322E794198C2F397DEF85F382FE9075A2A0E5F".to_string()));
        assert!(r < *CURVE_ORDER);
        let r = random_mod_order(Some(&mut rng));
        assert_eq!(r, BigNum::from_hex("0000000000000000000000000000000002EB5988AC48026ABEAF0206276C9D1158B5A2BE12EAFF2097A9AD8D8CFFD64D".to_string()));

        for _ in 0..30 {
            let r1 = random_mod_order::<ThreadRng>(None);
            let r2 = random_mod_order::<ThreadRng>(None);
            assert_ne!(r1, r2);
        }
    }

    #[test]
    fn hash_tests() {
        let seed = b"11111111111111111111111111111111";
        let mut rng = ChaChaRng::from_seed(*seed);

        let mut r = [0u8; 48];
        rng.fill_bytes(&mut r);

        let h = hash_mod_order(&r);
        assert_eq!(h, BigNum::from_hex("000000000000000000000000000000006E9AAA41F32B7DDC4C0CB1150B5414EA52C5234C329BF8B10E10C67F418A83FA".to_string()));
        assert!(h < *CURVE_ORDER);
        let h1 = hash_mod_order(&r);
        assert_eq!(h, h1);
        rng.fill_bytes(&mut r);
        let h = hash_mod_order(&r);
        assert_eq!(h, BIG::from_hex("0000000000000000000000000000000005A15BCF8E435724C15C8581587B2A408834236982E0405F047CADB2060B489B".to_string()));
        assert!(h < *CURVE_ORDER);
    }
}
