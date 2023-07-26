use ark_std::vec::Vec;
use ark_ec::{PairingEngine, ProjectiveCurve/*, AffineCurve*/};
use ark_ec::msm::VariableBaseMSM;
use ark_ff::{PrimeField, UniformRand, One, Zero}; // Field};
use ark_std::rand::Rng;

use rayon::prelude::*;

#[cfg(feature = "parallel")]

pub struct BasiliskSrs<E: PairingEngine> {
    pub srs: BasSrs<E>,
    pub pi: BasPi<E>,
}

pub struct BasSrs<E: PairingEngine> {
    pub xk_g1: Vec<<E as PairingEngine>::G1Projective>,
    pub x_g2: <E as PairingEngine>::G2Projective,
    pub u: <E as PairingEngine>::Fr,
}

pub struct BasPi<E: PairingEngine> {
    pub agg: BasPiAgg<E>,
    pub ind: BasPiInd<E>,
}

pub struct BasPiAgg<E: PairingEngine> {
    pub x_g1: <E as PairingEngine>::G1Projective,
}

pub struct BasPiInd<E: PairingEngine> {
    pub x_g1: <E as PairingEngine>::G1Projective,
    pub x_g2: <E as PairingEngine>::G2Projective,
}

#[inline]
/// Generates a random common reference string for
/// a circuit of size n.
pub fn basilisk_srs_gen<E,R>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, rng: &mut R, n: usize) -> BasiliskSrs<E> 
where
    E: PairingEngine,
    R: Rng,
{
    let n = n + 4;

    let x0 = E::Fr::rand(rng);
    let u0 = E::Fr::rand(rng); // doesn't have to be random

    let x0_g2 = g2_generator.mul(&x0.into_repr());
    let mut x0k_g1: Vec<E::G1Projective> = vec![g1_generator];
    let mut temp = x0;
    for _ in 1..n+1 {
        x0k_g1.push(g1_generator.mul(&temp.into_repr()));
        temp = temp * x0;
    }

    let pi_ind0 = BasPiInd {x_g1: x0k_g1[1], 
                            x_g2: x0_g2};
    let srs0 = BasSrs { xk_g1: x0k_g1, 
                        x_g2: x0_g2, 
                        u: u0};
    let pi_agg0 = BasPiAgg {x_g1: g1_generator.mul(&x0.into_repr())};
    let pi0 = BasPi {agg: pi_agg0, ind: pi_ind0};
    let crs = BasiliskSrs {srs: srs0, pi: pi0};

    return crs
}

#[inline]
pub fn basilisk_srs_update<E,R>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, rng: &mut R, srs: &BasSrs<E>, n: usize) -> BasiliskSrs<E>
where
    E: PairingEngine,
    R: Rng,
{
    let n = n+4;

    let x_bar = E::Fr::rand(rng);
    let u_bar = E::Fr::rand(rng); // doesn't need to be random

    let ui      = srs.u * u_bar;
    let xi_g2   = srs.x_g2.mul(&x_bar.into_repr());
    let mut xik_g1: Vec<E::G1Projective> = vec![g1_generator];
    let mut x_bar_pow = x_bar;
    for i in 1..n+1 { 
        xik_g1.push(srs.xk_g1[i].mul(&x_bar_pow.into_repr())); 
        x_bar_pow = x_bar_pow * x_bar;
    }

    let pi_indi = BasPiInd {x_g1: g1_generator.mul(&x_bar.into_repr()), 
                            x_g2: g2_generator.mul(&x_bar.into_repr())};
    let pi_aggi = BasPiAgg {x_g1: xik_g1[1]};
    let pii = BasPi {agg: pi_aggi, ind: pi_indi};
    let srsi = BasSrs {    xk_g1: xik_g1, 
                        x_g2: xi_g2, 
                        u: ui};
    let crs = BasiliskSrs {srs: srsi, pi: pii};
    return crs
    
}

/// uses multithreading to further improve performance
#[inline]
pub fn basilisk_srs_update_parallel<E,R>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, rng: &mut R, srs: &BasSrs<E>, n: usize) -> BasiliskSrs<E>
where
    E: PairingEngine,
    R: Rng,
{
    let n = n+4;

    let x_bar = E::Fr::rand(rng);
    let u_bar = E::Fr::rand(rng); // doesn't need to be random

    let ui      = srs.u * u_bar;
    let xi_g2   = srs.x_g2.mul(&x_bar.into_repr());
    
    let mut xk_xb: Vec<(E::G1Projective, E::Fr)> = vec![(srs.xk_g1[0], x_bar)];
    let mut xb = x_bar;
    for i in 1..n+1 {
        xk_xb.push((srs.xk_g1[i] , xb));
        xb *= x_bar;
    }

    let mut xik_g1 = xk_xb.into_par_iter().map(|xx| xx.0.mul(xx.1.into_repr())).collect::<Vec<_>>();
    xik_g1[0]= g1_generator;

    let pi_indi = BasPiInd {x_g1: g1_generator.mul(&x_bar.into_repr()),
                            x_g2: g2_generator.mul(&x_bar.into_repr())};
    let pi_aggi = BasPiAgg {x_g1: xik_g1[1]};
    let pii = BasPi {agg: pi_aggi, ind: pi_indi};
    let srsi = BasSrs {    xk_g1: xik_g1,
                        x_g2: xi_g2,
                        u: ui};
    let crs = BasiliskSrs {srs: srsi, pi: pii};
    return crs

}


pub fn basilisk_srs_verify<E:PairingEngine>(srs: &BasSrs<E>, pis: &Vec<BasPi<E>>, party: &str, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize) -> bool {

    match party {
        "prover"    => return verify_prover(srs, pis, g2_generator, n+4),
        "verifier"  => return verify_verifier(srs, pis, g1_generator, g2_generator, n+4),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_verifier<E:PairingEngine>(srs: &BasSrs<E>, pis: &Vec<BasPi<E>>, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize) -> bool {

    if pis.len() == 1 {             // for i = 0
        return true;
    }
    else {
        //check 1
        if pis[0].agg.x_g1 != pis[0].ind.x_g1 || srs.u == E::Fr::zero()  {
            return false
        }
        // check 2
        for i in 0..pis.len() {
            let eq2 = E::miller_loop([
                    (pis[i].ind.x_g1.into_affine().into(), g2_generator.into_affine().into()),
                    (g1_generator.into_affine().into(), (-pis[i].ind.x_g2.into_affine()).into())
                    ].iter());
            let test2 = E::final_exponentiation(&eq2).unwrap();
            if !(test2.is_one())    {return false;}
        }
        // check 3
        for i in 1..pis.len() {
            let eq3 = E::miller_loop([
                (pis[i].agg.x_g1.into_affine().into(), (-g2_generator.into_affine()).into()),
                (pis[i-1].agg.x_g1.into_affine().into(), pis[i].ind.x_g2.into_affine().into())
                ].iter());
            let test3 = E::final_exponentiation(&eq3).unwrap();
            if !(test3.is_one())    {return false;}
        }
        // check 4
        for k in 1..n+1 {
            let eq4 = E::miller_loop([
                (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
                (srs.xk_g1[k-1].into_affine().into(), srs.x_g2.into_affine().into())
                ].iter());
            let test4 = E::final_exponentiation(&eq4).unwrap();
            if !(test4.is_one())    {return false;}
        }
        return true
    }
}

fn verify_prover<E:PairingEngine>(srs: &BasSrs<E>, _pis: &Vec<BasPi<E>>, g2_generator: E::G2Projective, n: usize) -> bool {
    for k in 1..n+1 {
        let eq4 = E::miller_loop([
            (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
            (srs.xk_g1[k-1].into_affine().into(), srs.x_g2.into_affine().into())
            ].iter());
        let test4 = E::final_exponentiation(&eq4).unwrap();
        if !(test4.is_one())    {return false;}
    }
    if srs.u == E::Fr::zero() {return false;}
    return true
}

pub fn basilisk_srs_verify_batched<E:PairingEngine>(srs: &BasSrs<E>, pis: &Vec<BasPi<E>>, party: &str, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize, rand_st: &Vec< Vec<E::Fr> >, rand_h: &Vec<E::Fr>) -> bool {
    
    match party {
        "prover"    => return verify_prover_batched(srs, pis, g1_generator, g2_generator, n+4, rand_st,rand_h),
        "verifier"  => return verify_verifier_batched(srs, pis, g1_generator, g2_generator, n+4, rand_st, rand_h),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_prover_batched<E:PairingEngine>(srs: &BasSrs<E>, pis: &Vec<BasPi<E>>, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize, _rand_st: &Vec< Vec<E::Fr> >, rand_h: &Vec<E::Fr>) -> bool {
    
    let rand_alt: &[<E::Fr as PrimeField>::BigInt] = &rand_h.iter().map(|h| h.into_repr()).collect::<Vec<_>>()[..];
    let mut xk_1: Vec<E::G1Projective> = Vec::new();
    for j in 0..srs.xk_g1.len() {xk_1.push(srs.xk_g1[j]);}
    let xkg1: &[E::G1Affine] = &xk_1.into_par_iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let temp_l: E::G1Projective = pis[pis.len()-1].agg.x_g1 + VariableBaseMSM::multi_scalar_mul(&xkg1[2..n+1], rand_alt);
    let temp_r: E::G1Projective = g1_generator + VariableBaseMSM::multi_scalar_mul(&xkg1[1..n], rand_alt);
    
    let eq1 = E::miller_loop([
                (temp_l.into_affine().into(), (-g2_generator.into_affine()).into()),
                (temp_r.into_affine().into(), srs.x_g2.into_affine().into())
                ].iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();

    if !(test1.is_one() && srs.u != E::Fr::zero()) {
        return false;
    }
    return true
}


fn verify_verifier_batched<E:PairingEngine>(srs: &BasSrs<E>, pis: &Vec<BasPi<E>>, g1_generator: E::G1Projective, g2_generator: E::G2Projective, _n: usize, rand_st: &Vec<Vec<E::Fr>>, rand_h: &Vec<E::Fr>) -> bool {

    if pis.len() < 2 {return true}

    if !(pis[0].agg.x_g1 == pis[0].ind.x_g1) {return false}

    let mut xk_1: Vec<E::G1Projective> = Vec::new();
    for j in 0..srs.xk_g1.len() {xk_1.push(srs.xk_g1[j]);}

    let h: &[<E::Fr as PrimeField>::BigInt] = &rand_h.into_par_iter().map(|h| h.into_repr()).collect::<Vec<_>>()[..];
    let s: &[<E::Fr as PrimeField>::BigInt] = &rand_st.into_par_iter().map(|st| st[0].into()).collect::<Vec<_>>()[..];
    let t: &[<E::Fr as PrimeField>::BigInt] = &rand_st.into_par_iter().map(|st| st[1].into()).collect::<Vec<_>>()[..];

    let xbar_l: &[E::G1Affine] = &pis[1..].into_par_iter().map(|pi| pi.ind.x_g1.into_affine()).collect::<Vec<_>>()[..];
    let x_l: &[E::G1Affine] = &pis[1..].into_par_iter().map(|pi| pi.agg.x_g1.into_affine()).collect::<Vec<_>>()[..];

    let xkg1: &[E::G1Affine] = &xk_1.into_par_iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let xbar_r: &[E::G2Affine] = &pis[1..].into_par_iter().map(|pi| pi.ind.x_g2.into_affine()).collect::<Vec<_>>()[..];

    let temp_l: E::G1Projective = pis[0].ind.x_g1 + VariableBaseMSM::multi_scalar_mul(&xkg1[1..], h)
            + VariableBaseMSM::multi_scalar_mul(xbar_l, t) + VariableBaseMSM::multi_scalar_mul(x_l, s);

    let temp_r1: E::G2Projective = pis[0].ind.x_g2 + VariableBaseMSM::multi_scalar_mul(xbar_r, t);
    let temp_r3: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xkg1, h);

    let mut terms: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(pis.len()+3);  
    terms.push((temp_l.into_affine().into(),        (-g2_generator.into_affine()).into()));
    terms.push((g1_generator.into_affine().into(),  temp_r1.into_affine().into()));
    terms.push((temp_r3.into_affine().into(),       srs.x_g2.into_affine().into()));

    for j in 1..pis.len() {
        terms.push((pis[j-1].agg.x_g1.mul(s[j-1]).into_affine().into(), pis[j].ind.x_g2.into_affine().into() ));
    }

    let eq1 = E::miller_loop(terms.iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();
    if !(test1.is_one()) {return false}

    return true
} 

