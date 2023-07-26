//use crate::{
//    r1cs_to_qap::{LibsnarkReduction, R1CStoQAP},ProvingKey,Vec,VerifyingKey,
//};
use ark_std::vec::Vec;
use ark_ec::{PairingEngine, ProjectiveCurve/*, AffineCurve*/};
use ark_ff::{PrimeField, UniformRand, Field, One};// Field, Zero};
use ark_std::rand::Rng;
use ark_ec::msm::VariableBaseMSM;

#[cfg(feature = "parallel")]
//use rayon::prelude::*;

pub struct SonicSrs<E: PairingEngine> {
    pub srs: SonSrs<E>,
    pub pi: SonPi<E>,
}
pub struct SonSrs<E: PairingEngine> {
    pub xk_g1:  Vec<<E as PairingEngine>::G1Projective>,
    pub xk_g2:  Vec<<E as PairingEngine>::G2Projective>,
    pub axk_g1: Vec<<E as PairingEngine>::G1Projective>,
    pub axk_g2: Vec<<E as PairingEngine>::G2Projective>,
    //pub a_gT:   <E as PairingEngine>::Fr, // redundant
}
pub struct SonPi<E: PairingEngine> {
    pub agg: SonPiAgg<E>,
    pub ind: SonPiInd<E>,
}
pub struct SonPiAgg<E: PairingEngine> {
    pub x_g1:  <E as PairingEngine>::G1Projective,
//    pub x_g2:  <E as PairingEngine>::G2Projective,
    pub ax_g1: <E as PairingEngine>::G1Projective,
//    pub ax_g2: <E as PairingEngine>::G2Projective,
    pub a_g2:  <E as PairingEngine>::G2Projective,
}
pub struct SonPiInd<E: PairingEngine> {
    pub x_g1:  <E as PairingEngine>::G1Projective,
    pub x_g2:  <E as PairingEngine>::G2Projective,
    pub ax_g1: <E as PairingEngine>::G1Projective,
    pub ax_g2: <E as PairingEngine>::G2Projective,
    pub a_g2:  <E as PairingEngine>::G2Projective,
}

#[inline]
/// Generates a random common reference string for
/// a circuit.
pub fn sonic_srs_gen<E,R>(rng: &mut R, n: usize) -> SonicSrs<E> 
where
    E: PairingEngine,
    R: Rng,
{
    let g1_generator = E::G1Projective::rand(rng);
    let g2_generator = E::G2Projective::rand(rng);

    let x0 = E::Fr::rand(rng);
    let x0_inv = x0.inverse().unwrap();
    let a0 = E::Fr::rand(rng);


    let mut x0k_g1  = vec![g1_generator; 2*n+1];
    let mut x0k_g2 = vec![g2_generator; 2*n+1];
    let mut ax0k_g1 = vec![g1_generator; 2*n+1];
    let mut ax0k_g2 = vec![g2_generator.mul(a0.into_repr()); 2*n+1];
    for i in 1..n+1 {
        x0k_g1[n+i] = x0k_g1[n+i-1].mul(&x0.into_repr());
        x0k_g2[n+i] = x0k_g2[n+i-1].mul(&x0.into_repr());
        ax0k_g1[n+i] = x0k_g1[n+i].mul(&a0.into_repr());
        ax0k_g2[n+i] = x0k_g2[n+i].mul(&a0.into_repr());

        x0k_g1[n-i] = x0k_g1[n-i+1].mul(&x0_inv.into_repr());
        x0k_g2[n-i] = x0k_g2[n-i+1].mul(&x0_inv.into_repr());
        ax0k_g1[n-i] = x0k_g1[n-i].mul(&a0.into_repr());
        ax0k_g2[n-i] = x0k_g2[n-i].mul(&a0.into_repr());
    }
    
    let pia0 = SonPiAgg{x_g1:   x0k_g1[n+1], 
//                        x_g2:   x0k_g2[n+1],
                        ax_g1:  ax0k_g1[n+1],
//                        ax_g2:  x0k_g2[n+1].mul(a0.into_repr()),
                        a_g2:   ax0k_g2[n]  };
    let pii0 = SonPiInd{x_g1:   x0k_g1[n+1], 
                        x_g2:   x0k_g2[n+1],
                        ax_g1:  ax0k_g1[n+1],
                        ax_g2:  x0k_g2[n+1].mul(a0.into_repr()),
                        a_g2:   ax0k_g2[n]  };
    let pi0 = SonPi {   agg: pia0,
                        ind: pii0   };
    let srs0 = SonSrs {xk_g1:  x0k_g1, 
                        xk_g2:  x0k_g2, 
                        axk_g1:  ax0k_g1,
                        axk_g2:  ax0k_g2    };
    let srs = SonicSrs {srs: srs0, pi: pi0};

    return srs
}

#[inline]
/// Update the SRS for the new party i
pub fn sonic_srs_update<E,R>(rng: &mut R, srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, n: usize) -> SonicSrs<E>
where
    E: PairingEngine,
    R: Rng,
{
    let g1_generator = srs.xk_g1[n];
    let g2_generator = srs.xk_g2[n];
    
    let x_bar = E::Fr::rand(rng);
    let x_bar_inv = x_bar.inverse().unwrap();
    let a_bar = E::Fr::rand(rng);
    
    let mut xik_g1  = vec![g1_generator; 2*n+1];
    let mut xik_g2 = vec![g2_generator; 2*n+1];
    let mut axik_g1 = vec![g1_generator; 2*n+1];
    let mut axik_g2 = vec![srs.axk_g2[n].mul(&a_bar.into_repr()); 2*n+1];
  
    let mut x_bar_pow = x_bar;
    let mut x_bar_ipow = x_bar_inv;
    for i in 1..n+1 { 
        xik_g1[n+i] = srs.xk_g1[n+i].mul(&x_bar_pow.into_repr());
        xik_g2[n+i] = srs.xk_g2[n+i].mul(&x_bar_pow.into_repr());
        axik_g1[n+i] = srs.axk_g1[n+i].mul(&(x_bar_pow*a_bar).into_repr());
        axik_g2[n+i] = srs.axk_g2[n+i].mul(&(x_bar_pow*a_bar).into_repr());

        xik_g1[n-i] = srs.xk_g1[n-i].mul(&x_bar_ipow.into_repr());
        xik_g2[n-i] = srs.xk_g2[n-i].mul(&x_bar_ipow.into_repr());
        axik_g1[n-i] = srs.axk_g1[n-i].mul(&(x_bar_ipow*a_bar).into_repr());
        axik_g2[n-i] = srs.axk_g2[n-i].mul(&(x_bar_ipow*a_bar).into_repr());
        
        x_bar_pow = x_bar_pow * x_bar;
        x_bar_ipow = x_bar_ipow * x_bar_inv;

    }

    let pia = SonPiAgg{ x_g1:   xik_g1[n+1],
//                        x_g2:   xik_g2[n+1],
                        ax_g1:  axik_g1[n+1],
//                        ax_g2:  axik_g2[n+1],
                        a_g2:   pis[pis.len()-1].agg.a_g2.mul(&a_bar.into_repr())};
    let pii = SonPiInd{ x_g1:   g1_generator.mul(&x_bar.into_repr()), 
                        x_g2:   g2_generator.mul(&x_bar.into_repr()),
                        ax_g1:  g1_generator.mul(&(x_bar*a_bar).into_repr()),
                        ax_g2:  g2_generator.mul(&(x_bar*a_bar).into_repr()), 
                        a_g2:   g2_generator.mul(&a_bar.into_repr())  };
    let pi = SonPi {    agg: pia,
                        ind: pii};
    let srsi = SonSrs { xk_g1:  xik_g1, 
                        xk_g2:  xik_g2, 
                        axk_g1: axik_g1,
                        axk_g2: axik_g2  };
    let srs = SonicSrs {srs: srsi, pi: pi};
    return srs
    
}

pub fn sonic_srs_verify<E:PairingEngine>(srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, party: &str, n: usize) -> bool {
    match party {
        "prover"    => return verify_prover(srs, pis, n),
        "verifier"  => return verify_verifier(srs, pis, n),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_prover<E:PairingEngine>(srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, n: usize) -> bool {

    let g1_generator = srs.xk_g1[n];
    let g2_generator = srs.xk_g2[n];
    let i = pis.len();

    // for i = 0
    for k in 0..2*n+1 {         
        // CHECK A
        let eq1 = E::miller_loop([
            (srs.xk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
            (g1_generator.into_affine().into(), srs.xk_g2[k].into_affine().into()),
                    ].iter());
        let test1 = E::final_exponentiation(&eq1).unwrap();
        
        if !(test1.is_one()) {return false}

        if !(k==n || k==0) {
            // CHECK C
            let eq2 = E::miller_loop([
                (srs.axk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), srs.axk_g2[k].into_affine().into()),
            ].iter());
            let test2 = E::final_exponentiation(&eq2).unwrap();
            
            let eq3 = E::miller_loop([
                (srs.axk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
                (srs.xk_g1[k].into_affine().into(), pis[i-1].agg.a_g2.into_affine().into()),
                ].iter());
            let test3 = E::final_exponentiation(&eq3).unwrap();
            
            if !(test2.is_one() && test3.is_one()) {return false}
        }
        if !(k==0) {
            // CHECK B
            let eq4 = E::miller_loop([
                (srs.xk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
                (srs.xk_g1[k-1].into_affine().into(), srs.xk_g2[n+1].into_affine().into()),                ].iter());
            let test4 = E::final_exponentiation(&eq4).unwrap();
            
            if !(test4.is_one()) {return false}
        }
    }
    return true
}

fn verify_verifier<E:PairingEngine>(srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, n: usize) -> bool {

    let g1_generator = srs.xk_g1[n];
    let g2_generator = srs.xk_g2[n];
    let i = pis.len();

    if i < 2 {return true;}
    
    // check 1
    if !(pis[0].agg.x_g1 == pis[0].ind.x_g1 && /* pis[0].agg.x_g2 == pis[0].ind.x_g2 &&*/ pis[0].agg.ax_g1 == pis[0].ind.ax_g1 &&/* pis[0].agg.ax_g2 == pis[0].ind.ax_g2 &&*/ pis[0].agg.a_g2 == pis[0].ind.a_g2 ) {return false;}
    // check 2 
    for j in 0..i {
        // CHECK A
        let eq1 = E::miller_loop([
            (pis[j].ind.x_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
            (g1_generator.into_affine().into(), pis[j].ind.x_g2.into_affine().into()),
                    ].iter());
        let test1 = E::final_exponentiation(&eq1).unwrap();
        if !(test1.is_one()) {return false}

        // CHECK B
        let eq2 = E::miller_loop([
            (pis[j].ind.ax_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
            (pis[j].ind.x_g1.into_affine().into(), pis[j].ind.a_g2.into_affine().into()),
            ].iter());
        let test2 = E::final_exponentiation(&eq2).unwrap();
        let eq3 = E::miller_loop([
            (pis[j].ind.ax_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
            (g1_generator.into_affine().into(), pis[j].ind.ax_g2.into_affine().into()),
            ].iter());
        let test3 = E::final_exponentiation(&eq3).unwrap();
        if !(test2.is_one() && test3.is_one()) {return false}
    }
    for j in 1..i {
        // CHECK C
        let eq4a = E::miller_loop([
            (pis[j].agg.x_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
            (pis[j-1].agg.x_g1.into_affine().into(), pis[j].ind.x_g2.into_affine().into()),
            ].iter());
        let test4a = E::final_exponentiation(&eq4a).unwrap();
//        let eq4b = E::miller_loop([
//            (pis[j].agg.x_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
//            (g1_generator.into_affine().into(), pis[j].agg.x_g2.into_affine().into()),
//            ].iter());
//        let test4b = E::final_exponentiation(&eq4b).unwrap();
        
        if !(test4a.is_one()/* && test4b.is_one()*/) {return false}

        // CHECK D
//        let eq5a = E::miller_loop([
//            (pis[j].agg.ax_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
//            (g1_generator.into_affine().into(), pis[j].agg.ax_g2.into_affine().into()),
//            ].iter());
//        let test5a = E::final_exponentiation(&eq5a).unwrap();
        let eq5b = E::miller_loop([
            (pis[j].agg.ax_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
            (pis[j].agg.x_g1.into_affine().into(), pis[j].agg.a_g2.into_affine().into()),
            ].iter());
        let test5b = E::final_exponentiation(&eq5b).unwrap();
        let eq5c = E::miller_loop([
            (pis[j].agg.ax_g1.into_affine().into(),(-g2_generator.into_affine()).into()),
            (pis[j-1].agg.ax_g1.into_affine().into(), pis[j].ind.ax_g2.into_affine().into()),
            ].iter());
        let test5c = E::final_exponentiation(&eq5c).unwrap();
        if !(/*test5a.is_one() &&*/ test5b.is_one() && test5c.is_one()) {return false}
    }

    for k in 0..2*n+1 {
        //check E
        let eq6 = E::miller_loop([
            (srs.xk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
            (g1_generator.into_affine().into(), srs.xk_g2[k].into_affine().into()),
            ].iter());
        let test6 = E::final_exponentiation(&eq6).unwrap();
        if !(test6.is_one()) {return false;}

        // CHECK H
        if !(k==0) {
            let eq7 = E::miller_loop([
                    (srs.xk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
                    (srs.xk_g1[k-1].into_affine().into(), srs.xk_g2[n+1].into_affine().into()),
                    ].iter());
            let test7 = E::final_exponentiation(&eq7).unwrap();
            if !(test7.is_one()) {return false;}
        }

        // CHECK I
        if !(k==n) {
            let eq8a = E::miller_loop([
                    (srs.axk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
                    (g1_generator.into_affine().into(), srs.axk_g2[k].into_affine().into()),
                    ].iter());
            let test8a = E::final_exponentiation(&eq8a).unwrap();
            
            let eq8b = E::miller_loop([
                    (srs.axk_g1[k].into_affine().into(),(-g2_generator.into_affine()).into()),
                    (srs.xk_g1[k].into_affine().into(), pis[i-1].agg.a_g2.into_affine().into()),
                    ].iter());
            let test8b = E::final_exponentiation(&eq8b).unwrap();
            if !(test8a.is_one() && test8b.is_one()) { return false}
        }
    }
    return true
}


pub fn sonic_srs_verify_batched<E:PairingEngine>(srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, party: &str, n: usize, rand_r: &Vec< Vec<E::Fr> >, rand_t: &Vec<Vec<E::Fr>>) -> bool {

    match party {
        "prover"    => return verify_prover_batched(srs, pis, n, rand_r,rand_t),
        "verifier"  => return verify_verifier_batched(srs, pis, n, rand_r, rand_t),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

pub fn verify_prover_batched<E:PairingEngine>(srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, n: usize, _rand_r: &Vec<Vec<E::Fr>>, rand_t: &Vec<Vec<E::Fr>>) -> bool {
    let g1_generator = srs.xk_g1[n];
    let g2_generator = srs.xk_g2[n];
    let i = pis.len();
   
    let t: &[<E::Fr as PrimeField>::BigInt] = &rand_t[0].iter().map(|t| t.into_repr()).collect::<Vec<_>>()[..];
    let xkg1_l: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let xkg1_r: &[E::G2Affine] = &srs.xk_g2.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let temp_l: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xkg1_l, t);
    let temp_r: E::G2Projective = VariableBaseMSM::multi_scalar_mul(xkg1_r, t);

    let eq1 = E::miller_loop([
                (temp_l.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), temp_r.into_affine().into())
                ].iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();

    if !(test1.is_one()) {
        return false;}
    

    let temp_l: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&xkg1_l[1..], &t[1..]);
    let temp_r: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&xkg1_l[0..2*n], &t[1..]);

    let eq2 = E::miller_loop([
                (temp_l.into_affine().into(), (-g2_generator.into_affine()).into()),
                (temp_r.into_affine().into(), srs.xk_g2[n+1].into_affine().into())
                ].iter());
    let test2 = E::final_exponentiation(&eq2).unwrap();

    if !(test2.is_one()) {
        return false;}
    
    let t_hat: &[<E::Fr as PrimeField>::BigInt] = &rand_t[1].iter().map(|t| t.into_repr()).collect::<Vec<_>>()[..];
    //let mut t_h: Vec<<E::Fr as PrimeField>::BigInt> = t_hat[0..n].to_vec();
    //t_h.append(&mut t_hat[n+1..].to_vec());
    //let t_hat : &[<E::Fr as PrimeField>::BigInt] = &t_h[..];

    let axkg1: &[E::G1Affine] = &srs.axk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let mut ax: Vec<E::G1Affine> = axkg1[0..n].to_vec();
    ax.append(&mut axkg1[n+1..].to_vec());
    let axkg1: &[E::G1Affine] = &ax[..];

    let axkg2: &[E::G2Affine] = &srs.axk_g2.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let mut ax: Vec<E::G2Affine> = axkg2[0..n].to_vec();
    ax.append(&mut axkg2[n+1..].to_vec());
    let axkg2: &[E::G2Affine] = &ax[..];
    
    let xkg1: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let mut x: Vec<E::G1Affine> = xkg1[0..n].to_vec();
    x.append(&mut xkg1[n+1..].to_vec());
    let xkg1: &[E::G1Affine] = &x[..];
    
    let temp_a: E::G1Projective = VariableBaseMSM::multi_scalar_mul(axkg1, t_hat);
    let temp_b: E::G2Projective = VariableBaseMSM::multi_scalar_mul(axkg2, t_hat);
    let temp_c: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xkg1,  t_hat);

    let eq3a = E::miller_loop([
                (temp_a.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), temp_b.into_affine().into())
                ].iter());
    let test3a = E::final_exponentiation(&eq3a).unwrap();
    let eq3b = E::miller_loop([
                (temp_a.into_affine().into(), (-g2_generator.into_affine()).into()),
                (temp_c.into_affine().into(), pis[i-1].agg.a_g2.into_affine().into())
                ].iter());
    let test3b = E::final_exponentiation(&eq3b).unwrap();

    if !(test3a.is_one() && test3b.is_one()) {
        return false;
    }
    return true;
}

pub fn verify_verifier_batched<E:PairingEngine>(srs: &SonSrs<E>, pis: &Vec<SonPi<E>>, n: usize, rand_r: &Vec<Vec<E::Fr>>, rand_t: &Vec<Vec<E::Fr>>) -> bool {
    let g1_generator = srs.xk_g1[n];
    let g2_generator = srs.xk_g2[n];
    let i = pis.len();

    if i < 2 {return true;}
    // check 1
    if !(pis[0].agg.x_g1 == pis[0].ind.x_g1 /*&& pis[0].agg.x_g2 == pis[0].ind.x_g2*/ && pis[0].agg.ax_g1 == pis[0].ind.ax_g1 /*&& pis[0].agg.ax_g2 == pis[0].ind.ax_g2*/ && pis[0].agg.a_g2 == pis[0].ind.a_g2 ) {return false;}
    
    // check 2 
    let r1: &[<E::Fr as PrimeField>::BigInt] = &rand_r[0].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let xbar1: &[E::G1Affine] = &pis.iter().map(|pi| pi.ind.x_g1.into_affine()).collect::<Vec<_>>()[..];
    let xbar2: &[E::G2Affine] = &pis.iter().map(|pi| pi.ind.x_g2.into_affine()).collect::<Vec<_>>()[..];
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xbar1, r1);
    let tempb: E::G2Projective = VariableBaseMSM::multi_scalar_mul(xbar2, r1);

    let eq1 = E::miller_loop([
                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), tempb.into_affine().into())
                ].iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();

    if !(test1.is_one()) {return false;}

    // check 3 
    let r2: &[<E::Fr as PrimeField>::BigInt] = &rand_r[1].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let axbar1: &[E::G1Affine] = &pis.iter().map(|pi| pi.ind.ax_g1.into_affine()).collect::<Vec<_>>()[..];
    let axbar2: &[E::G2Affine] = &pis.iter().map(|pi| pi.ind.ax_g2.into_affine()).collect::<Vec<_>>()[..];
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(axbar1, r2);
    let tempb: E::G2Projective = VariableBaseMSM::multi_scalar_mul(axbar2, r2);

    let mut terms2b: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i+1);
    terms2b.push((tempa.into_affine().into(), (-g2_generator.into_affine()).into()));
    for j in 0..i {
        terms2b.push((pis[j].ind.x_g1.mul(r2[j]).into_affine().into(), pis[j].ind.a_g2.into_affine().into() ));
    }

    let eq2a = E::miller_loop([
                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), tempb.into_affine().into())
                ].iter());
    let test2a = E::final_exponentiation(&eq2a).unwrap();
    let eq2b = E::miller_loop(terms2b.iter());
    let test2b = E::final_exponentiation(&eq2b).unwrap();

    if !(test2a.is_one() && test2b.is_one()) {return false;}

    // check 4
    let r3: &[<E::Fr as PrimeField>::BigInt] = &rand_r[2].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let x1: &[E::G1Affine] = &pis.iter().map(|pi| pi.agg.x_g1.into_affine()).collect::<Vec<_>>()[..];
//    let x2: &[E::G2Affine] = &pis.iter().map(|pi| pi.agg.x_g2.into_affine()).collect::<Vec<_>>()[..];
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&x1[1..], &r3[1..]);
//    let tempb: E::G2Projective = VariableBaseMSM::multi_scalar_mul(&x2[1..], &r3[1..]);

    let mut terms3b: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i+1);
    terms3b.push((tempa.into_affine().into(), (-g2_generator.into_affine()).into()));
    for j in 1..i {
        terms3b.push((pis[j-1].agg.x_g1.mul(r3[j]).into_affine().into(), pis[j].ind.x_g2.into_affine().into()));
    }

//    let eq3a = E::miller_loop([
//                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
//                (g1_generator.into_affine().into(), tempb.into_affine().into())
//                ].iter());
//    let test3a = E::final_exponentiation(&eq3a).unwrap();
    let eq3b = E::miller_loop(terms3b.iter());
    let test3b = E::final_exponentiation(&eq3b).unwrap();

    if !(/*test3a.is_one() && */test3b.is_one()) {return false;}

    // check 5
    let r4: &[<E::Fr as PrimeField>::BigInt] = &rand_r[3].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let ax1: &[E::G1Affine] = &pis.iter().map(|pi| pi.agg.ax_g1.into_affine()).collect::<Vec<_>>()[..];
//    let ax2: &[E::G2Affine] = &pis.iter().map(|pi| pi.agg.ax_g2.into_affine()).collect::<Vec<_>>()[..];
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&ax1[1..], &r4[1..]);
//    let tempb: E::G2Projective = VariableBaseMSM::multi_scalar_mul(&ax2[1..], &r4[1..]);

    let mut terms4b: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i+1);
    terms4b.push((tempa.into_affine().into(), (-g2_generator.into_affine()).into()));
    for j in 1..i {
        terms4b.push((pis[j].agg.x_g1.mul(r4[j]).into_affine().into(), pis[j].agg.a_g2.into_affine().into()));
    }
    let mut terms4c: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i+1);
    terms4c.push((tempa.into_affine().into(), (-g2_generator.into_affine()).into()));
    for j in 1..i {
        terms4c.push((pis[j-1].agg.ax_g1.mul(r4[j]).into_affine().into(), pis[j].ind.ax_g2.into_affine().into()));
    }


//    let eq4a = E::miller_loop([
//                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
//                (g1_generator.into_affine().into(), tempb.into_affine().into())
//                ].iter());
//    let test4a = E::final_exponentiation(&eq4a).unwrap();
    let eq4b = E::miller_loop(terms4b.iter());
    let test4b = E::final_exponentiation(&eq4b).unwrap();
    let eq4c = E::miller_loop(terms4c.iter());
    let test4c = E::final_exponentiation(&eq4c).unwrap();

    if !(/*test4a.is_one() &&*/ test4b.is_one() && test4c.is_one()) {return false;}

    // check 6 
    let t1: &[<E::Fr as PrimeField>::BigInt] = &rand_t[0].iter().map(|t| t.into_repr()).collect::<Vec<_>>()[..];
    let xk1: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let xk2: &[E::G2Affine] = &srs.xk_g2.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xk1, t1);
    let tempb: E::G2Projective = VariableBaseMSM::multi_scalar_mul(xk2, t1);

    let eq5 = E::miller_loop([
                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), tempb.into_affine().into())
                ].iter());
    let test5 = E::final_exponentiation(&eq5).unwrap();
    
    if !(test5.is_one()) {return false;}

    // check 7 
    let xk1: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&xk1[1..], &t1[1..]);
    let tempb: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&xk1[0..], &t1[1..]);

    let eq6 = E::miller_loop([
                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
                (tempb.into_affine().into(), srs.xk_g2[n+1].into_affine().into())
                ].iter());
    let test6 = E::final_exponentiation(&eq6).unwrap();

    if !(test6.is_one()) {return false;}

    // check 8
    let t2: &[<E::Fr as PrimeField>::BigInt] = &rand_t[1].iter().map(|t| t.into_repr()).collect::<Vec<_>>()[..];
    let mut t22: Vec<<E::Fr as PrimeField>::BigInt> = t2[0..n].to_vec();
    t22.append(&mut t2[n+1..].to_vec());
    let t2: &[<E::Fr as PrimeField>::BigInt] = &t22[..];
    
    let axk1: &[E::G1Affine] = &srs.axk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let mut ax: Vec<E::G1Affine> = axk1[0..n].to_vec();
    ax.append(&mut axk1[n+1..].to_vec());
    let axk1: &[E::G1Affine] = &ax[..];
    
    let axk2: &[E::G2Affine] = &srs.axk_g2.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let mut ax: Vec<E::G2Affine> = axk2[0..n].to_vec();
    ax.append(&mut axk2[n+1..].to_vec());
    let axk2: &[E::G2Affine] = &ax[..];
    
    //let xk1: &[E::G1Affine] = &xk1[0..n].append(xk1[n+1..]);
    let mut x: Vec<E::G1Affine> = xk1[0..n].to_vec();
    x.append(&mut xk1[n+1..].to_vec());
    let xk1: &[E::G1Affine] = &x[..];
    
    let tempa: E::G1Projective = VariableBaseMSM::multi_scalar_mul(axk1, t2);
    let tempb: E::G2Projective = VariableBaseMSM::multi_scalar_mul(axk2, t2);
    let tempc: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xk1, t2);


    let eq7a = E::miller_loop([
                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), tempb.into_affine().into())
                ].iter());
    let test7a = E::final_exponentiation(&eq7a).unwrap();
    let eq7b = E::miller_loop([
                (tempa.into_affine().into(), (-g2_generator.into_affine()).into()),
                (tempc.into_affine().into(), pis[i-1].agg.a_g2.into_affine().into())
                ].iter());
    let test7b = E::final_exponentiation(&eq7b).unwrap();

    if !(test7a.is_one() && test7b.is_one()) {return false;}

    return true;
}
