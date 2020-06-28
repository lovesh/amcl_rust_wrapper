extern crate rand;
extern crate sha3;

use rand::{CryptoRng, RngCore};

use crate::constants::{
    BarrettRedc_FM_k, BarrettRedc_FM_u, BarrettRedc_FM_v, BigNumBits, CurveOrder, Ell,
    FieldElement_SIZE, FieldModulus, DNLEN,
};
use crate::types::{BigNum, DoubleBigNum, FP};
use amcl::hmac;
use amcl::rand::RAND;

use sha3::digest::{ExtendableOutput, Input, XofReader};
use sha3::Shake256;

/// Hash message and return output of size equal to curve modulus. Uses SHAKE to hash the message.
pub fn hash_msg(msg: &[u8]) -> [u8; FieldElement_SIZE] {
    let mut hasher = Shake256::default();
    hasher.input(&msg);
    let mut h: [u8; FieldElement_SIZE] = [0; FieldElement_SIZE];
    hasher.xof_result().read(&mut h);
    h
}

pub fn get_seeded_RNG_with_rng<R: RngCore + CryptoRng>(entropy_size: usize, rng: &mut R) -> RAND {
    // initialize from at least 128 byte string of raw random entropy
    let mut entropy = vec![0; entropy_size];
    rng.fill_bytes(&mut entropy.as_mut_slice());
    get_RAND(entropy_size, entropy.as_slice())
}

pub fn get_seeded_RNG(entropy_size: usize) -> RAND {
    let mut entropy = vec![0; entropy_size];
    let mut rng = rand::thread_rng();
    rng.fill_bytes(&mut entropy.as_mut_slice());
    get_RAND(entropy_size, entropy.as_slice())
}

fn get_RAND(entropy_size: usize, entropy: &[u8]) -> RAND {
    let mut r = RAND::new();
    r.clean();
    r.seed(entropy_size, &entropy);
    r
}

/// Hash to one or more field elements, largely copied from miracl/core.
/// Used in hash to curve
pub fn hash_to_field(
    hash: usize,
    hlen: usize,
    domain_separation_tag: &[u8],
    msg: &[u8],
    output: &mut [FP],
    out_count: usize,
) {
    let mut okm: [u8; 256] = [0; 256];
    let mut fd: [u8; 128] = [0; 128];

    hmac::xmd_expand(
        hash,
        hlen,
        &mut okm,
        *Ell * out_count,
        &domain_separation_tag,
        &msg,
    );
    for i in 0..out_count {
        for j in 0..*Ell {
            fd[j] = okm[*Ell * i + j];
        }
        // Use Barrett reduction since the FieldModulus is fixed for a curve
        let fp = FP::new_big(&barrett_reduction(
            &DoubleBigNum::frombytes(&fd[0..*Ell]),
            &FieldModulus,
            *BarrettRedc_FM_k,
            &*BarrettRedc_FM_u,
            &*BarrettRedc_FM_v,
        ));
        output[i] = fp;
    }
}

// set x = x mod 2^m
fn mod2m(x: &mut DoubleBigNum, m: usize) {
    // the limb to update
    let limb_idx = m / BigNumBits;
    // number of bits in limb at `limb_index` to keep the same as in `x`
    let bt = m % BigNumBits;
    // mask is all 1s
    let msk = (1 << bt) - 1;
    x.w[limb_idx] &= msk;
    // Set the remaining bits to 0
    for i in limb_idx + 1..DNLEN {
        x.w[i] = 0
    }
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
    let mut q1 = DoubleBigNum::new_copy(x);
    q1.shr(k - 1);
    // Above right shift will convert q from DBIG to BIG
    let q1 = BigNum::new_dcopy(&q1);

    let q2 = BigNum::mul(&q1, &u);

    // q3 = floor(q2 / 2^{k+1})
    let mut q3 = DoubleBigNum::new_copy(&q2);
    q3.shr(k + 1);
    let q3 = BigNum::new_dcopy(&q3);

    // r1 = x % 2^{k+1}
    let mut r1 = DoubleBigNum::new_copy(x);
    mod2m(&mut r1, k + 1);
    let r1 = BigNum::new_dcopy(&r1);

    // r2 = (q3 * modulus) % 2^{k+1}
    let mut r2 = BigNum::mul(&q3, modulus);
    mod2m(&mut r2, k + 1);
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
    r.norm();
    r
}

// Reducing BigNum for comparison with `rmod`
fn __barrett_reduction__(x: &BigNum, modulus: &BigNum, k: usize, u: &BigNum, v: &BigNum) -> BigNum {
    // q1 = floor(x / 2^{k-1})
    let mut q1 = x.clone();
    q1.shr(k - 1);

    let q2 = BigNum::mul(&q1, &u);

    // q3 = floor(q2 / 2^{k+1})
    let mut q3 = DoubleBigNum::new_copy(&q2);
    q3.shr(k + 1);
    let q3 = BigNum::new_dcopy(&q3);

    // r1 = x % 2^{k+1}
    let mut r1 = x.clone();
    r1.mod2m(k + 1);

    // r2 = (q3 * modulus) % 2^{k+1}
    let mut r2 = BigNum::mul(&q3, modulus);
    mod2m(&mut r2, k + 1);
    let r2 = BigNum::new_dcopy(&r2);

    // if r1 > r2, r = r1 - r2 else r = r1 - r2 + v
    // Since negative numbers are not supported, use r2 - r1. This holds since r = r1 - r2 + v = v - (r2 - r1)
    let diff = BigNum::comp(&r1, &r2);
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
    r.norm();
    r
}

/// For a modulus returns
/// k = number of bits in modulus
/// u = floor(2^2k / modulus)
/// v = 2^(k+1)
pub fn barrett_reduction_params(modulus: &BigNum) -> (usize, BigNum, BigNum) {
    let k = modulus.nbits();

    // u = floor(2^2k/CurveOrder)
    let mut u = DoubleBigNum::new();
    u.w[0] = 1;
    // `u.shl(2*k)` crashes, so perform shl(k) twice
    u.shl(k);
    u.shl(k);
    // div returns floored value
    let u = u.div(&CurveOrder);

    // v = 2^(k+1)
    let mut v = BigNum::new_int(1isize);
    v.shl(k + 1);

    (k, u, v)
}

/// Takes a collection of collections and pads them with the given element such that all collections are of same length.
/// Returns the new length. Mutates the collections.
#[macro_export]
macro_rules! pad_collection {
    ( $coll:expr, $pad:expr ) => {{
        if $coll.len() > 0 {
            // Find maximum length of a collection
            let mut max_length = $coll[0].len();
            for i in 1..$coll.len() {
                if max_length < $coll[i].len() {
                    max_length = $coll[i].len();
                }
            }
            for i in 0..$coll.len() {
                let l = $coll[i].len();
                if l < max_length {
                    // append might not unroll and hence perform slower
                    //$coll[i].append(&mut vec![$pad; max_length - l]);
                    for _ in 0..max_length - l {
                        $coll[i].push($pad);
                    }
                }
            }

            max_length
        } else {
            0
        }
    }};
}

/// Checks the given string begins and ends with the given characters, i.e. "bounded" and removes
/// the beginning and ending characters. Throws error if not "bounded"
#[macro_export]
macro_rules! unbound_bounded_string {
    ( $string:ident, $start_bound:expr, $end_bound:expr, $error:expr ) => {
        // Need $string as "<$start_bound>x,y<$end_bound>", so if $start_bound was '(' and
        // $end_bound was ')', string should be (foo..bar....)
        if !$string.starts_with($start_bound)
            || !$string.ends_with($end_bound)
            || !$string.len() < 2
        {
            return Err($error);
        }
        // Remove bounds
        $string.remove(0);
        $string.pop();
    };
}

/// Expects the string to be a concatenation of 2 strings of equal length separated by comma and
/// returns tuple of 2 strings
#[macro_export]
macro_rules! split_string_to_2_tuple {
    ( $string:ident, $error:expr ) => {{
        // Not using split to avoid cloning strings
        let mut str_2 = $string.split_off($string.len() / 2);
        if !str_2.starts_with(',') {
            return Err($error);
        }
        // Remove ','
        str_2.remove(0);
        ($string, str_2)
    }};
}

/// Expects the string to be a concatenation of 3 strings of equal length separated by commas and
/// returns tuple of 3 strings
#[macro_export]
macro_rules! split_string_to_3_tuple {
    ( $string:ident, $error:expr ) => {{
        // Not using split to avoid cloning strings
        let mut bc = $string.split_off($string.len() / 3);
        if !bc.starts_with(',') {
            return Err($error);
        }
        bc.remove(0);

        let mut c = bc.split_off(bc.len() / 2);
        if !c.starts_with(',') {
            return Err($error);
        }
        c.remove(0);
        ($string, bc, c)
    }};
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
    use std::time::{Duration, Instant};

    #[test]
    fn timing_fp_big() {
        // TODO: Compare adding raw BIGs and FieldElement to check the overhead of the abstraction
        let count = 100;
        let elems: Vec<_> = (0..count).map(|_| FieldElement::random()).collect();
        let bigs: Vec<_> = elems.iter().map(|f| f.to_bignum()).collect();
        let fs: Vec<_> = bigs.iter().map(|b| FP::new_big(&b)).collect();
        let mut res_mul = BIG::new_int(1 as isize);
        let mut start = Instant::now();
        for b in &bigs {
            res_mul = BigNum::modmul(&res_mul, &b, &CurveOrder);
        }
        println!(
            "Multiplication time for {} BIGs = {:?}",
            count,
            start.elapsed()
        );

        let mut res_mul = FP::new_int(1 as isize);
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
            i.invmodp(&CurveOrder);
            inverses_b.push(i);
        }
        println!("Inverse time for {} BIGs = {:?}", count, start.elapsed());
        for i in 0..count {
            let r = BigNum::modmul(&inverses_b[i], &bigs[i], &CurveOrder);
            assert_eq!(BigNum::comp(&r, &BigNum::new_int(1 as isize)), 0);
        }

        start = Instant::now();
        for f in &fs {
            let mut i = f.clone();
            i.inverse(None);
            inverses_f.push(i);
        }
        println!("Inverse time for {} FPs = {:?}", count, start.elapsed());
        for i in 0..count {
            let mut c = inverses_f[i].clone();
            c.mul(&fs[i]);
            assert!(c.equals(&FP::new_int(1 as isize)));
        }

        // Fixme: add in FP crashes while adding 100 elems
        let c = 50;
        start = Instant::now();
        let mut r = bigs[0];
        for i in 0..c {
            r.add(&bigs[i]);
            r.rmod(&CurveOrder);
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
        //let (k, u, v) = barrett_reduction_params(&CurveOrder);
        let (k, u, v) = (
            *constants::BarrettRedc_k,
            *constants::BarrettRedc_u,
            *constants::BarrettRedc_v,
        );
        let mut xs = vec![];
        let mut reduced1 = vec![];
        let mut reduced2 = vec![];
        let mut rng = rand::thread_rng();
        let count = 1000;
        for _ in 0..count {
            let a: u32 = rng.gen();
            let s = BigNum::new_int(a as isize);
            let _x = CurveOrder.minus(&s);
            xs.push(BigNum::mul(&_x, &_x));
        }

        let mut start = Instant::now();
        for x in &xs {
            let r = barrett_reduction(&x, &CurveOrder, k, &u, &v);
            reduced1.push(r);
        }
        println!("Barrett time = {:?}", start.elapsed());

        start = Instant::now();
        for x in &xs {
            let mut y = DoubleBigNum::new_copy(x);
            let z = y.dmod(&CurveOrder);
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
            *constants::BarrettRedc_k,
            *constants::BarrettRedc_u,
            *constants::BarrettRedc_v,
        );
        let count = 100;
        let elems: Vec<_> = (0..count).map(|_| FieldElement::random()).collect();
        let bigs: Vec<_> = elems.iter().map(|f| f.to_bignum()).collect();

        let mut sum = bigs[0].clone();
        let mut start = Instant::now();
        for i in 0..count {
            sum = BigNum::plus(&sum, &bigs[i]);
            sum.rmod(&CurveOrder)
        }
        println!("rmod time = {:?}", start.elapsed());

        let mut sum_b = bigs[0].clone();
        start = Instant::now();
        for i in 0..count {
            sum_b = BigNum::plus(&sum_b, &bigs[i]);
            sum_b = __barrett_reduction__(&sum_b, &CurveOrder, k, &u, &v)
        }
        println!("Barrett time = {:?}", start.elapsed());

        assert_eq!(BigNum::comp(&sum, &sum_b), 0)
    }
}
