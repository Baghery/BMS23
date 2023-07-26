#![warn(unused)]
#![deny(
    trivial_casts,
    trivial_numeric_casts,
    variant_size_differences,
    stable_features,
    non_shorthand_field_patterns,
    renamed_and_removed_lints,
    private_in_public,
    unsafe_code
)]
use ark_ec::PairingEngine;
use ark_ff::UniformRand;
use ark_std::test_rng;
use ark_std::rand::Rng;
use uu_snarks::*;

use ark_bls12_381::Bls12_381;
use ark_std::str::*;


fn get_random_bounded<E: PairingEngine>(bound: usize) -> E::Fr {
    // HARDCODED FOR bound = 80  
    
    let mut rand_gen = test_rng();
    let mut rand_str = rand_gen.gen_range(1..10).to_string().to_owned();
    for _ in 1..24 { // 24 = log10(2**80)
        rand_str = rand_str + &rand_gen.gen_range(0..10).to_string();
    }
    
    let rand = match E::Fr::from_str(&rand_str) {
        Ok(r)   => r,
        Err(_e)  => get_random_bounded::<E>(bound),
    };
    return rand
}

fn find_adversary<E: PairingEngine>(g1_generator: E::G1Projective, g2_generator: E::G2Projective, srss: &Vec<BasSrs<E>>, pis: &Vec<BasPi<E>>, n: usize, left: usize, right: usize) -> usize {
    if right > pis.len()  {return 666}
    let mid = (left+right)/2;

    let mut rand_st : Vec::<Vec::<E::Fr>> = Vec::new(); // ![vec![E::Fr::rand(rng)]];
    let mut rand_h : Vec::<E::Fr> = Vec::new();
    for _ in 0..n+5 {rand_h.push(get_random_bounded::<E>(80));}
    for _ in 0..mid {rand_st.push(vec!(get_random_bounded::<E>(80), get_random_bounded::<E>(80)));}

    if left+2 == right {
        let check = basilisk_srs_verify_batched(&srss[mid-1], &pis[..mid], "prover", g1_generator, g2_generator, n, &rand_st, &rand_h);
        if check {return left+1}
        else {return left}
    }

    
    let check_mid: bool = basilisk_srs_verify_batched(&srss[mid-1], &pis[..mid], "verifier", g1_generator, g2_generator, n, &rand_st, &rand_h);
    
    if check_mid {
        return find_adversary::<E>(g1_generator, g2_generator, &srss, &pis, n, mid, right);
    }
    else {
        return find_adversary::<E>(g1_generator, g2_generator, &srss, &pis, n, left, mid);
    }
}

fn identifiable_security<E: PairingEngine>(adversary:usize, n_updates: usize, c_size: usize)
{

    let rng = &mut test_rng();

    let g1_generator = E::G1Projective::rand(rng);
    let g2_generator = E::G2Projective::rand(rng);

    let mut pis: Vec<BasPi<E>> = Vec::new();
    let mut srss: Vec<BasSrs<E>> = Vec::new();

    let mut rand_st : Vec::<Vec::<E::Fr>> = Vec::new(); // ![vec![E::Fr::rand(rng)]];
    let mut rand_h : Vec::<E::Fr> = Vec::new();
    for _ in 0..c_size+5 {rand_h.push(get_random_bounded::<E>(80));}
    for _ in 1..n_updates{rand_st.push(vec!(get_random_bounded::<E>(80), get_random_bounded::<E>(80)));}

    let mut srsj = basilisk_srs_gen(g1_generator, g2_generator, rng, c_size);
    pis.push(srsj.pi);
    srss.push(srsj.srs);

    for j in 1..n_updates {

        if j == adversary {
            srsj = basilisk_srs_update_parallel(g1_generator, g2_generator, rng, &srss[srss.len()-1], c_size);
            pis.push(srsj.pi);
            srsj.srs.xk_g1[2] = srsj.srs.xk_g1[1];
            srss.push(srsj.srs);
        }
        else {
            srsj = basilisk_srs_update_parallel(g1_generator, g2_generator, rng, &srss[j-1], c_size);
            pis.push(srsj.pi);
            srss.push(srsj.srs);
        }
    }
    
    let check: bool = basilisk_srs_verify_batched(&srss[n_updates-1], &pis, "verifier", g1_generator, g2_generator, c_size, &rand_st, &rand_h);
    if !check {let eve = find_adversary::<E>(g1_generator, g2_generator, &srss,&pis,c_size,0,n_updates);
        assert!(eve == adversary);
        for j in eve..n_updates-1 {
            srsj = basilisk_srs_update_parallel(g1_generator, g2_generator, rng, &srss[j-1], c_size);
            pis[j] = srsj.pi;
            srss[j]= srsj.srs;
        }
        assert!(basilisk_srs_verify_batched(&srss[n_updates-2], &pis[..n_updates-1], "verifier", g1_generator, g2_generator, c_size, &rand_st, &rand_h));
    }
}

fn main() {
    let n_iters = 1;
    let n_updates = 8;
    let c_size = 10;
   
    for ad in 1..n_updates {
        identifiable_security::<Bls12_381>(ad, n_updates, c_size);
    }
}
