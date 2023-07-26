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

use num_bigint::{BigUint,RandBigInt};
use ark_ff::{BigInteger, PrimeField, BigInteger128};

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
        Err(e)  => get_random_bounded::<E>(bound),
    };
    return rand
}

fn basilisk_timing<E: PairingEngine>(_n_iters:usize, n_updates: usize, c_size: usize)
{
    let rng = &mut test_rng();
    let mut randgen = rand::thread_rng();

    let g1_generator = E::G1Projective::rand(rng);
    let g2_generator = E::G2Projective::rand(rng);

    let mut pis: Vec<BasPi<E>> = Vec::new();

    let mut rand_st : Vec::<Vec::<E::Fr>> = Vec::new(); // ![vec![E::Fr::rand(rng)]];
    let mut rand_h : Vec::<E::Fr> = Vec::new();
    for _ in 0..c_size+4 {
        rand_h.push(get_random_bounded::<E>(80));
    }
    for _ in 0..n_updates {
        rand_st.push(vec!(get_random_bounded::<E>(80), get_random_bounded::<E>(80)));
    }

    let time_sg = ark_std::time::Instant::now();
    let mut srsj = basilisk_srs_gen(g1_generator, g2_generator, rng, c_size);
    print!("basilisk_gen: \t\t{}\ni: \tsu\t\tsvp\t\tsvv\n", time_sg.elapsed().as_nanos());
    pis.push(srsj.pi);

    for i in 1..n_updates {
        let time_su = ark_std::time::Instant::now();
        srsj = basilisk_srs_update_parallel(g1_generator, g2_generator, rng, &srsj.srs, c_size);
        let su = time_su.elapsed().as_nanos();
        pis.push(srsj.pi);

        let time_svb_p = ark_std::time::Instant::now();
        assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "prover", g1_generator, g2_generator,c_size, &rand_st, &rand_h));
        let svb_p = time_svb_p.elapsed().as_nanos();
        let time_svb_v = ark_std::time::Instant::now();
        assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "verifier", g1_generator, g2_generator, c_size, &rand_st, &rand_h));
        let svb_v = time_svb_v.elapsed().as_nanos();
        
        print!("{}:\t{}\t{}\t{}\n",i,su,svb_p,svb_v);
    
    
    }

    // regular timing 
    for i in 0..n_iters {
        let time_sv_p = ark_std::time::Instant::now();
        assert!(basilisk_srs_verify(&srsj.srs, &pis, "prover", g1_generator, g2_generator, c_size));
        print!("basilisk_verification (prover) {}: \t{}\n",i, time_sv_p.elapsed().as_nanos());
        
        let time_sv_v = ark_std::time::Instant::now();
        assert!(basilisk_srs_verify(&srsj.srs, &pis, "verifier", g1_generator, g2_generator, c_size));
        print!("basilisk_verification (verifier) {}: \t{}\n",i, time_sv_v.elapsed().as_nanos());
    }

    // batched timing
    let mut rand_st : Vec::<Vec::<E::Fr>> = Vec::new(); // ![vec![E::Fr::rand(rng)]];
    let mut rand_h : Vec::<E::Fr> = Vec::new();
    for _ in 0..c_size {
        rand_h.push(E::Fr::rand(rng));
    }
    for _ in 0..n_updates {
        rand_st.push(vec!(E::Fr::rand(rng), E::Fr::rand(rng)));
    }

    let time_svb_p = ark_std::time::Instant::now();
    assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "prover", g1_generator, g2_generator,c_size, &rand_st, &rand_h));
    print!("basilisk_verification (prover): \t{}\n", time_svb_p.elapsed().as_nanos());
    
    let time_svb_v = ark_std::time::Instant::now();
    assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "verifier", g1_generator, g2_generator, c_size, &rand_st, &rand_h));
    print!("basilisk_verification (verifier): \t{}\n", time_svb_v.elapsed().as_nanos());

}

fn main() {
    let n_iters = 1;
    let n_updates = 6;
    let c_size = 5000;
    
    print!("\ni = {}, N = {}\n", n_updates, c_size);
    basilisk_timing::<Bls12_381>(n_iters, n_updates, c_size);
}
