use ark_std::vec::Vec;
use ark_ec::{PairingEngine, ProjectiveCurve/*, AffineCurve*/};
use ark_ec::msm::VariableBaseMSM;
use ark_ff::{PrimeField, UniformRand, One};// Field, Zero};
use ark_std::rand::Rng;

#[cfg(feature = "parallel")]

pub struct MarlinSrs<E: PairingEngine> {
    pub srs: MarSrs<E>,
    pub pi: MarPi<E>,
}

pub struct MarSrs<E: PairingEngine> {
    pub xk_g1:  Vec<<E as PairingEngine>::G1Projective>,
    pub gxk_g1: Vec<<E as PairingEngine>::G1Projective>,
    pub x_g2:   <E as PairingEngine>::G2Projective,
    pub gx_g2:  <E as PairingEngine>::G2Projective,
}

pub struct MarPi<E: PairingEngine> {
    pub agg: MarPiAgg<E>,
    pub ind: MarPiInd<E>,
}

pub struct MarPiAgg<E: PairingEngine> {
    pub g_g1: <E as PairingEngine>::G1Projective,
    pub gx_g1: <E as PairingEngine>::G1Projective,
    pub x_g2: <E as PairingEngine>::G2Projective,
}

pub struct MarPiInd<E: PairingEngine> {
    pub x_g1: <E as PairingEngine>::G1Projective,
    pub g_g1: <E as PairingEngine>::G1Projective,
    pub x_g2: <E as PairingEngine>::G2Projective,
    pub gx_g2: <E as PairingEngine>::G2Projective,
}

#[inline]
/// Generates a random common reference string for
/// a circuit.
pub fn marlin_srs_gen<E,R>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, rng: &mut R, n: usize) -> MarlinSrs<E> 
where
    E: PairingEngine,
    R: Rng,
{
    let x0 = E::Fr::rand(rng);
    let g0 = E::Fr::rand(rng);


    let mut xk_g1: Vec<E::G1Projective> = vec![g1_generator];
    let mut gxk_g1: Vec<E::G1Projective> = vec![g1_generator.mul(g0.into_repr())];
    let mut xk = x0;
    for _ in 1..n+1 {
        xk_g1.push(xk_g1[0].mul(xk.into_repr()));
        gxk_g1.push(gxk_g1[0].mul(xk.into_repr()));
        xk *= x0;
    }

    let x_g2 = g2_generator.mul(x0.into_repr());
    let gx_g2 = g2_generator.mul((g0*x0).into_repr());

    let pi_agg = MarPiAgg { g_g1:   gxk_g1[0],
                            gx_g1:  gxk_g1[1],
                            x_g2:   x_g2};
    let pi_ind = MarPiInd { x_g1:   xk_g1[1],
                            g_g1:   gxk_g1[0],
                            x_g2:   x_g2,
                            gx_g2:  gx_g2};
    let srs = MarSrs      { xk_g1:  xk_g1,
                            gxk_g1: gxk_g1,
                            x_g2:   x_g2, 
                            gx_g2:  gx_g2};
    let pi = MarPi        {agg: pi_agg, ind: pi_ind};
    let crs = MarlinSrs    {srs: srs, pi: pi};

    return crs
}

#[inline]
pub fn marlin_srs_update<E,R>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, rng: &mut R, srs: &MarSrs<E>, pis: &Vec<MarPi<E>>, n: usize) -> MarlinSrs<E>
where
    E: PairingEngine,
    R: Rng,
{
    let x_bar = E::Fr::rand(rng);
    let g_bar = E::Fr::rand(rng);
    let i = pis.len();

    let mut xk_g1: Vec<E::G1Projective> = vec![g1_generator];
    let mut gxk_g1: Vec<E::G1Projective> = vec![pis[i-1].agg.g_g1.mul(g_bar.into_repr())];
    let mut xk = x_bar;
    let mut gxk = x_bar*g_bar;
    for k in 1..n+1 {
        xk_g1.push(srs.xk_g1[k].mul(xk.into_repr()));
        gxk_g1.push(srs.gxk_g1[k].mul(gxk.into_repr()));
        xk *= x_bar;
        gxk *= x_bar;
    }
    let x_g2 = srs.x_g2.mul(x_bar.into_repr());
    let gx_g2 = srs.gx_g2.mul((g_bar*x_bar).into_repr());


    let pi_agg = MarPiAgg { g_g1:   gxk_g1[0],
                            gx_g1:  gxk_g1[1],
                            x_g2:   x_g2};
    let pi_ind = MarPiInd { x_g1:   g1_generator.mul(x_bar.into_repr()),
                            g_g1:   g1_generator.mul(g_bar.into_repr()),
                            x_g2:   g2_generator.mul(x_bar.into_repr()),
                            gx_g2:  g2_generator.mul((g_bar*x_bar).into_repr())};
    let srs = MarSrs       {xk_g1:  xk_g1,
                            gxk_g1: gxk_g1,
                            x_g2:   x_g2,
                            gx_g2:  gx_g2};
    let pi = MarPi        {agg: pi_agg, ind: pi_ind};
    let crs = MarlinSrs    {srs: srs, pi: pi};
    
    return crs
    
}

pub fn marlin_srs_verify<E:PairingEngine>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, srs: &MarSrs<E>, pis: &Vec<MarPi<E>>, party: &str, n: usize) -> bool {

    match party {
        "prover"    => return verify_prover(srs, g2_generator, n),
        "verifier"  => return verify_verifier(srs, pis, g1_generator, g2_generator, n),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_prover<E:PairingEngine>(srs: &MarSrs<E>, g2_generator: E::G2Projective, n: usize) -> bool {
    for k in 1..n+1 {
        let eq1 = E::miller_loop([
            (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
            (srs.xk_g1[k-1].into_affine().into(), srs.x_g2.into_affine().into())
            ].iter());
        let test1 = E::final_exponentiation(&eq1).unwrap();
        if !(test1.is_one())    {return false;}
            
        let eq2 = E::miller_loop([
            (srs.gxk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
            (srs.gxk_g1[k-1].into_affine().into(), srs.x_g2.into_affine().into())
            ].iter());
        let test2 = E::final_exponentiation(&eq2).unwrap();
        if !(test2.is_one())    {return false;}
        
    }
    return true
}

fn verify_verifier<E:PairingEngine>(srs: &MarSrs<E>, pis: &Vec<MarPi<E>>, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize) -> bool {
    
    let i = pis.len();
    if i < 2 {             // for i = 0
        return true;
    }
    else {
        //check 1
        if !(pis[0].agg.g_g1 == pis[0].ind.g_g1 && pis[0].agg.x_g2 == pis[0].ind.x_g2)  {
            return false
        }
        // check 2
        for i in 0..pis.len() {
            let eq2a = E::miller_loop([
                    (pis[i].ind.x_g1.into_affine().into(), g2_generator.into_affine().into()),
                    (g1_generator.into_affine().into(), (-pis[i].ind.x_g2.into_affine()).into())
                    ].iter());
            let test2a = E::final_exponentiation(&eq2a).unwrap();
            let eq2b = E::miller_loop([
                    (pis[i].ind.g_g1.into_affine().into(), pis[i].ind.x_g2.into_affine().into()),
                    (g1_generator.into_affine().into(), (-pis[i].ind.gx_g2.into_affine()).into())
                    ].iter());
            let test2b = E::final_exponentiation(&eq2b).unwrap();
            if !(test2a.is_one() && test2b.is_one())    {
                return false;}
        }
        // check 3
        for i in 1..pis.len() {
            let eq3a = E::miller_loop([
                ((-g1_generator.into_affine()).into(), pis[i].agg.x_g2.into_affine().into()),
                (pis[i].ind.x_g1.into_affine().into(), pis[i-1].agg.x_g2.into_affine().into())
                ].iter());
            let test3a = E::final_exponentiation(&eq3a).unwrap();
            let eq3b = E::miller_loop([
                (pis[i].agg.gx_g1.into_affine().into(), (-g2_generator.into_affine()).into()),
                (pis[i-1].agg.gx_g1.into_affine().into(), pis[i].ind.gx_g2.into_affine().into())
                ].iter());
            let test3b = E::final_exponentiation(&eq3b).unwrap();
            let eq3c = E::miller_loop([
                (pis[i].agg.gx_g1.into_affine().into(), (-g2_generator.into_affine()).into()),
                (pis[i].agg.g_g1.into_affine().into(), pis[i].agg.x_g2.into_affine().into())
                ].iter());
            let test3c = E::final_exponentiation(&eq3c).unwrap();
            if !(test3a.is_one() && test3b.is_one() && test3c.is_one())    {
                return false;}
        }
        // check 4
        for k in 1..n+1 {
            let eq4 = E::miller_loop([
                (srs.xk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
                (srs.xk_g1[k-1].into_affine().into(), srs.x_g2.into_affine().into())
                ].iter());
            let test4 = E::final_exponentiation(&eq4).unwrap();
            if !(test4.is_one())    {
                return false;}
            
            let eq5 = E::miller_loop([
                (srs.gxk_g1[k].into_affine().into(), (-g2_generator.into_affine()).into()),
                (srs.gxk_g1[k-1].into_affine().into(), srs.x_g2.into_affine().into())
                ].iter());
            let test5 = E::final_exponentiation(&eq5).unwrap();
            if !(test5.is_one())    {
                return false;}
        }

        let eq6 = E::miller_loop([
            (srs.gxk_g1[1].into_affine().into(), (-g2_generator.into_affine()).into()),
            (g1_generator.into_affine().into(), srs.gx_g2.into_affine().into())
            ].iter());
        let test6 = E::final_exponentiation(&eq6).unwrap();
        if !(test6.is_one())    {
            return false;}

        return true
    }
}

pub fn marlin_srs_verify_batched<E:PairingEngine>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, srs: &MarSrs<E>, pis: &Vec<MarPi<E>>, party: &str, n: usize, rand_r: &Vec< Vec<E::Fr> >, rand_t: &Vec<E::Fr>) -> bool {
    
    match party {
        "prover"    => return verify_prover_batched(srs, pis, g1_generator, g2_generator, n, rand_t),
        "verifier"  => return verify_verifier_batched(srs, pis, g1_generator, g2_generator, n, rand_r, rand_t),
        _           => panic!("party should be either \"prover\" or \"verifier\"."),
    }
}

fn verify_prover_batched<E:PairingEngine>(srs: &MarSrs<E>, _pis: &Vec<MarPi<E>>, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize, rand_t: &Vec<E::Fr>) -> bool {
    
    let t: &[<E::Fr as PrimeField>::BigInt] = &rand_t.iter().map(|t| t.into_repr()).collect::<Vec<_>>()[..];
    let xkg1: &[E::G1Affine] = &srs.xk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    let gxkg1: &[E::G1Affine] = &srs.gxk_g1.iter().map(|x| x.into_affine()).collect::<Vec<_>>()[..];
    
    let mut base_l = xkg1[1..].to_vec();
    base_l.append(&mut gxkg1[1..].to_vec());
    let mut base_r = xkg1[..n].to_vec();
    base_r.append(&mut gxkg1[..n].to_vec());
    let msm_l: E::G1Projective = srs.gxk_g1[1] + VariableBaseMSM::multi_scalar_mul(&base_l[..], &t[..]);
    let msm_r: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&base_r[..], &t[..]);
    
    let eq1 = E::miller_loop([
                (msm_l.into_affine().into(), (-g2_generator.into_affine()).into()),
                (msm_r.into_affine().into(), srs.x_g2.into_affine().into()),
                (g1_generator.into_affine().into(), srs.gx_g2.into_affine().into())
                ].iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();

    if !(test1.is_one()) {
        return false;
    }
    return true
}

fn verify_verifier_batched<E:PairingEngine>(srs: &MarSrs<E>, pis: &Vec<MarPi<E>>, g1_generator: E::G1Projective, g2_generator: E::G2Projective, n: usize, rand_r: &Vec<Vec<E::Fr>>, rand_t: &Vec<E::Fr>) -> bool {

    let i = pis.len();
    if i < 2 {return true}

    if !(pis[0].agg.g_g1 == pis[0].ind.g_g1 && pis[0].agg.x_g2 == pis[0].ind.x_g2) {
        print!("fail 0\n");
        return false}

    let r1: &[<E::Fr as PrimeField>::BigInt] = &rand_r[0].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let r2: &[<E::Fr as PrimeField>::BigInt] = &rand_r[1].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let r3: &[<E::Fr as PrimeField>::BigInt] = &rand_r[2].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let r4: &[<E::Fr as PrimeField>::BigInt] = &rand_r[3].iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];
    let t: &[<E::Fr as PrimeField>::BigInt] = &rand_t.iter().map(|r| r.into_repr()).collect::<Vec<_>>()[..];

    let xbj_1: &[E::G1Affine] = &pis.iter().map(|pi| pi.ind.x_g1.into_affine()).collect::<Vec<_>>()[..];
    let xbj_2: &[E::G2Affine] = &pis.iter().map(|pi| pi.ind.x_g2.into_affine()).collect::<Vec<_>>()[..];
    let gxbj_2: &[E::G2Affine] = &pis.iter().map(|pi| pi.ind.gx_g2.into_affine()).collect::<Vec<_>>()[..];
    let xj_2: &[E::G2Affine] = &pis.iter().map(|pi| pi.agg.x_g2.into_affine()).collect::<Vec<_>>()[..];
    let gxj_1: &[E::G1Affine] = &pis.iter().map(|pi| pi.agg.gx_g1.into_affine()).collect::<Vec<_>>()[..];

    let xk_1: &[E::G1Affine] = &srs.xk_g1.iter().map(|xk| xk.into_affine()).collect::<Vec<_>>()[..];
    let gxk_1: &[E::G1Affine] = &srs.gxk_g1.iter().map(|gxk| gxk.into_affine()).collect::<Vec<_>>()[..];
    
    let msm1_a: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&xbj_1[..], r1);
    let msm1_b: E::G2Projective = VariableBaseMSM::multi_scalar_mul(xbj_2, r1);
    let eq1 = E::miller_loop([
                (msm1_a.into_affine().into(), (-g2_generator.into_affine()).into()),
                (g1_generator.into_affine().into(), msm1_b.into_affine().into())
                ].iter());
    let test1 = E::final_exponentiation(&eq1).unwrap();
    if !(test1.is_one()) {
        print!("fail 1 {}\n", test1.is_one());
        return false}
    
    let mut terms2: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i+1);
    let msm2_a: E::G2Projective = VariableBaseMSM::multi_scalar_mul(gxbj_2,r2);
    terms2.push(( (-g1_generator.into_affine()).into(), msm2_a.into_affine().into() ));
    for j in 0..i {
        terms2.push(( pis[j].ind.g_g1.mul(r2[j]).into_affine().into(), 
                      pis[j].ind.x_g2.into_affine().into() ));
    }
    let eq2 = E::miller_loop(terms2.iter());
    let test2 = E::final_exponentiation(&eq2).unwrap();
    if !(test2.is_one()) {
        print!("fail 2 {}\n", test2.is_one());
        return false}

    let mut terms3: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i);
    let msm3_a: E::G2Projective = VariableBaseMSM::multi_scalar_mul(&xj_2[1..],&r3[1..]);
    terms3.push(( (-g1_generator.into_affine()).into(), msm3_a.into_affine().into() ));
    for j in 1..i {
        terms3.push(( pis[j].ind.x_g1.mul(r3[j]).into_affine().into(), 
                      pis[j-1].agg.x_g2.into_affine().into() ));
    }
    let eq3 = E::miller_loop(terms3.iter());
    let test3 = E::final_exponentiation(&eq3).unwrap();
    if !(test3.is_one()) {
        print!("fail 3 {}\n", test3.is_one());
        return false}

    let mut terms4a: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i);
    let mut terms4b: Vec<(E::G1Prepared, E::G2Prepared)> = Vec::with_capacity(i);
    let msm4_a: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&gxj_1[1..],&r4[1..]);
    terms4a.push(( msm4_a.into_affine().into(), (-g2_generator.into_affine()).into()));
    terms4b.push(( msm4_a.into_affine().into(), (-g2_generator.into_affine()).into()));
    for j in 1..i {
        terms4a.push((pis[j-1].agg.gx_g1.mul(r4[j]).into_affine().into(), 
                      pis[j].ind.gx_g2.into_affine().into() ));
        terms4b.push((pis[j].agg.g_g1.mul(r4[j]).into_affine().into(), 
                      pis[j].agg.x_g2.into_affine().into() ));
    }
    let eq4a = E::miller_loop(terms4a.iter());
    let eq4b = E::miller_loop(terms4b.iter());
    let test4a = E::final_exponentiation(&eq4a).unwrap();
    let test4b = E::final_exponentiation(&eq4b).unwrap();
    if !(test4a.is_one() && test4b.is_one()) {
        print!("fail 4 {} {}\n", test4a.is_one(), test4b.is_one());
        return false}

    let mut xkgxk_1 = xk_1[1..].to_vec();
    xkgxk_1.append(&mut gxk_1[1..].to_vec());
    let mut xkgxk0_1 = xk_1[..n].to_vec();
    xkgxk0_1.append(&mut gxk_1[..n].to_vec());
    let msm5_a: E::G1Projective = srs.gxk_g1[1] +  VariableBaseMSM::multi_scalar_mul(&xkgxk_1[..],t);
    let msm5_b: E::G1Projective = VariableBaseMSM::multi_scalar_mul(&xkgxk0_1[..], t);
    let eq5 = E::miller_loop([
        (msm5_a.into_affine().into(),  (-g2_generator.into_affine()).into()),
        (msm5_b.into_affine().into(), srs.x_g2.into_affine().into()),
        (g1_generator.into_affine().into(), srs.gx_g2.into_affine().into()), 
        ].iter());
    let test5 = E::final_exponentiation(&eq5).unwrap();
    if !(test5.is_one()) {
        print!("fail 5 {}\n", test5.is_one());
        return false}
    
    return true
} 

