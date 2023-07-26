use ark_std::vec::Vec;
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ec::msm::VariableBaseMSM;
use ark_ff::{PrimeField, UniformRand, One};// Field, Zero};
use ark_std::rand::Rng;

#[cfg(feature = "parallel")]

pub struct LunarLiteSrs<E: PairingEngine> {
    pub srs: LunSrs<E>,
    pub pi: LunPi<E>,
}

pub struct LunSrs<E: PairingEngine> {
    pub xk_g1:  Vec<<E as PairingEngine>::G1Projective>,
    pub xk_g2:  Vec<<E as PairingEngine>::G2Projective>,
}

pub struct LunPi<E: PairingEngine> {
    pub agg: LunPiAgg<E>,
    pub ind: LunPiInd<E>,
}

pub struct LunPiAgg<E: PairingEngine> {
    pub x_g1:  <E as PairingEngine>::G1Projective,
}

pub struct LunPiInd<E: PairingEngine> {
    pub x_g1:  <E as PairingEngine>::G1Projective,
    pub x_g2:  <E as PairingEngine>::G2Projective,
}

#[inline]
pub fn lunarlite_srs_gen<E,R>(rng: &mut R, n: usize) -> LunarLiteSrs<E> 
where
    E: PairingEngine,
    R: Rng,
{
    let g1_generator = E::G1Projective::rand(rng);
    let g2_generator = E::G2Projective::rand(rng);

    let x0 = E::Fr::rand(rng);

    let mut x0k_g1  = vec![g1_generator; n+1];
    let mut x0k_g2 = vec![g2_generator; n+1];
    for i in 1..n+1 {
        x0k_g1[i] = x0k_g1[i-1].mul(&x0.into_repr());
        x0k_g2[i] = x0k_g2[i-1].mul(&x0.into_repr());
    }

    let pi_ind0 = LunPiInd {x_g1:   x0k_g1[1], 
                            x_g2:   x0k_g2[1],};
    let pi_agg0 = LunPiAgg {x_g1:   x0k_g1[1]};
    let pi0 = LunPi {agg: pi_agg0, ind: pi_ind0};
    let srs0 = LunSrs {xk_g1:  x0k_g1, 
                        xk_g2:  x0k_g2};
    let srs = LunarLiteSrs {srs: srs0, pi: pi0};

    return srs
}

#[inline]
/// Update the SRS for the new party i
pub fn lunarlite_srs_update<E,R>(rng: &mut R, srs: &LunSrs<E>, _pis: &Vec<LunPi<E>>, n: usize) -> LunarLiteSrs<E>
where
    E: PairingEngine,
    R: Rng,
{
    let g1_generator = srs.xk_g1[0];
    let g2_generator = srs.xk_g2[0];
    
    let x_bar = E::Fr::rand(rng);
    
    let mut xik_g1  = vec![g1_generator; n+1];
    let mut xik_g2 = vec![g2_generator; n+1];
  
    let mut x_bar_pow = x_bar;
    for i in 1..n+1 { 
        xik_g1[i] = srs.xk_g1[i].mul(&x_bar_pow.into_repr());
        xik_g2[i] = srs.xk_g2[i].mul(&x_bar_pow.into_repr());

        x_bar_pow = x_bar_pow * x_bar;
    }

    let pia = LunPiAgg {x_g1:   xik_g1[1]}; 
    let pii = LunPiInd {x_g1:   g1_generator.mul(&x_bar.into_repr()), 
                        x_g2:   g2_generator.mul(&x_bar.into_repr())};
    let pi = LunPi {agg: pia, ind: pii};
    let srsi = LunSrs { xk_g1:  xik_g1, 
                        xk_g2:  xik_g2};
    let srs = LunarLiteSrs {srs: srsi, pi: pi};
    return srs
    
}

pub fn lunarlite_srs_verify<E:PairingEngine>(srs: &LunSrs<E>, pis: &Vec<LunPi<E>>, party: &str, n: usize) -> bool {

    match party { //true: Prover, false: Verifier
        "prover"    => return verify_prover(srs, pis, n),
        "verifier"  => return verify_verifier(srs, pis, n),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_prover<E:PairingEngine>(srs: &LunSrs<E>, pis: &Vec<LunPi<E>>, n: usize) -> bool {
    let g1_generator = srs.xk_g1[0];
    let g2_generator = srs.xk_g2[0];
    for k in 1..n+1 {
        let eq1 = E::miller_loop([
            (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
            (srs.xk_g1[k-1].into_affine().into(), srs.xk_g2[1].into_affine().into())
            ].iter());
        let eq2 = E::miller_loop([
            (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
            (g1_generator.into_affine().into(), srs.xk_g2[k].into_affine().into())
            ].iter());
        let test1 = E::final_exponentiation(&eq1).unwrap();
        let test2 = E::final_exponentiation(&eq2).unwrap();

        if !(test1.is_one() && test2.is_one())    {print!("false 0\n{},{}\n",test1.is_one(),test2.is_one());return false;}
    }
    return true
}

fn verify_verifier<E:PairingEngine>(srs: &LunSrs<E>, pis: &Vec<LunPi<E>>, n: usize) -> bool {

    let g1_generator = srs.xk_g1[0];
    let g2_generator = srs.xk_g2[0];
    
    if pis.len() == 1 {             // for i = 0
        return true;
    }
    else {
        //check 1
        if pis[0].agg.x_g1 != pis[0].ind.x_g1  {
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
                (srs.xk_g1[k-1].into_affine().into(), srs.xk_g2[1].into_affine().into())
                ].iter());
            let test4 = E::final_exponentiation(&eq4).unwrap();
            let eq5 = E::miller_loop([
                (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), srs.xk_g2[k].into_affine().into())
                ].iter());
            let test5 = E::final_exponentiation(&eq5).unwrap();
            if !(test4.is_one() && test5.is_one())    {return false;}
        }
        return true
    }
}

pub fn lunarlite_srs_verify_batched<E:PairingEngine>(srs: &LunSrs<E>, pis: &Vec<LunPi<E>>, party: &str, n: usize, rand_st: &Vec< Vec<E::Fr> >, rand_h: &Vec<E::Fr>) -> bool {

    match party {
        "prover"    => return verify_prover_batched(srs, pis, n, rand_st,rand_h),
        "verifier"  => return verify_verifier_batched(srs, pis, n, rand_st, rand_h),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_prover_batched<E:PairingEngine>(srs: &LunSrs<E>, pis: &Vec<LunPi<E>>, n: usize, _rand_st: &Vec< Vec<E::Fr> >, rand_h: &Vec<E::Fr>) -> bool {

    let g1_generator = srs.xk_g1[0];
    let g2_generator = srs.xk_g2[0];
    
    let rand_alt: &[<E::Fr as PrimeField>::BigInt] = &rand_h.iter().map(|h| h.into_repr()).collect::<Vec<_>>()[..];
    let xkg1_l: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let xkg1_m: &[E::G2Affine] = &srs.xk_g2.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];

    let temp_l: E::G1Projective = pis[pis.len()-1].agg.x_g1 + VariableBaseMSM::multi_scalar_mul(&xkg1_l[2..], rand_alt);
    let temp_r: E::G1Projective = g1_generator              + VariableBaseMSM::multi_scalar_mul(&xkg1_l[1..n], rand_alt);
    let temp_m: E::G2Projective = srs.xk_g2[1] + VariableBaseMSM::multi_scalar_mul(&xkg1_m[2..], rand_alt);

    let eq1 = E::miller_loop([
                (temp_l.into_affine().into(), (-g2_generator.into_affine()).into()),
                (temp_r.into_affine().into(), srs.xk_g2[1].into_affine().into())
                ].iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();
    let eq2 = E::miller_loop([
                (temp_l.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), temp_m.into_affine().into())
                ].iter());
    let test2 = E::final_exponentiation(&eq2).unwrap();

    if !(test1.is_one() && test2.is_one()) {
        print!("false vb\n{},{}\n",test1.is_one(),test2.is_one());
        return false;
    }
    return true
}

fn verify_verifier_batched<E:PairingEngine>(srs: &LunSrs<E>, pis: &Vec<LunPi<E>>, _n: usize, rand_st: &Vec<Vec<E::Fr>>, rand_h: &Vec<E::Fr>) -> bool {

    let g1_generator = srs.xk_g1[0];
    let g2_generator = srs.xk_g2[0];
    
    if pis.len() < 2 {return true}

    if !(pis[0].agg.x_g1 == pis[0].ind.x_g1) {print!("false 0\n");return false}

    let h: &[<E::Fr as PrimeField>::BigInt] = &rand_h.iter().map(|h| h.into_repr()).collect::<Vec<_>>()[..];
    let h2: &[<E::Fr as PrimeField>::BigInt] = &rand_h.iter().map(|h| (*h+h).into_repr()).collect::<Vec<_>>()[..];
    let s: &[<E::Fr as PrimeField>::BigInt] = &rand_st.iter().map(|st| st[0].into()).collect::<Vec<_>>()[..];
    let t: &[<E::Fr as PrimeField>::BigInt] = &rand_st.iter().map(|st| st[1].into()).collect::<Vec<_>>()[..];

    let xkg1_l: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let xbar_l: &[E::G1Affine] = &pis[1..].iter().map(|pi| pi.ind.x_g1.into_affine()).collect::<Vec<_>>()[..];
    let x_l: &[E::G1Affine] = &pis[1..].iter().map(|pi| pi.agg.x_g1.into_affine()).collect::<Vec<_>>()[..];

    let xkg2_r: &[E::G2Affine] = &srs.xk_g2.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let xbar_r: &[E::G2Affine] = &pis[1..].iter().map(|pi| pi.ind.x_g2.into_affine()).collect::<Vec<_>>()[..];


    let temp_l: E::G1Projective = pis[0].ind.x_g1 + VariableBaseMSM::multi_scalar_mul(&xkg1_l[1..], h2)
            + VariableBaseMSM::multi_scalar_mul(xbar_l, t) + VariableBaseMSM::multi_scalar_mul(x_l, s);

    let temp_r1: E::G2Projective = pis[0].ind.x_g2 + VariableBaseMSM::multi_scalar_mul(xbar_r, t)
            + VariableBaseMSM::multi_scalar_mul(&xkg2_r[1..], h);
    let temp_r3: E::G1Projective = VariableBaseMSM::multi_scalar_mul(xkg1_l, h);

    let mut terms: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(pis.len()+3);
    terms.push((temp_l.into_affine().into(),        (-g2_generator.into_affine()).into()));
    terms.push((g1_generator.into_affine().into(),  temp_r1.into_affine().into()));
    terms.push((temp_r3.into_affine().into(),       srs.xk_g2[1].into_affine().into()));

    for j in 1..pis.len() {
        terms.push((pis[j-1].agg.x_g1.mul(s[j-1]).into_affine().into(), pis[j].ind.x_g2.into_affine().into() ));
    }

    let eq1 = E::miller_loop(terms.iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();
    if !(test1.is_one()) {print!("false 1\n");return false}
    return true
}
