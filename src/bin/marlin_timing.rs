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
use uu_snarks::marlin_srs::*;
use ark_bls12_377::Bls12_377;

fn marlin_timing<E: PairingEngine>(n_iters:usize, n_updates: usize, c_size: usize)
{
    let rng = &mut test_rng();

    let g1_generator = E::G1Projective::rand(rng);
    let g2_generator = E::G2Projective::rand(rng);
  
    let mut pis: Vec<MarPi<E>> = Vec::new();
    let time_sg = ark_std::time::Instant::now();
    let mut srsj = marlin_srs_gen(g1_generator, g2_generator, rng, c_size);
    print!("marlin_gen: \t{}\n", time_sg.elapsed().as_nanos());
    pis.push(srsj.pi);

    for _ in 1..n_updates {
        let time_su = ark_std::time::Instant::now();
        srsj = marlin_srs_update(g1_generator, g2_generator, rng,&srsj.srs,&pis,c_size);
        print!("marlin_update: \t{}\n", time_su.elapsed().as_nanos());
        pis.push(srsj.pi);
    }
    
    // regular timing
    for i in 0..n_iters {
        let time_sv_p = ark_std::time::Instant::now();
        assert!(marlin_srs_verify(g1_generator, g2_generator, &srsj.srs, &pis, "prover", c_size));
        print!("marlin_verification (prover) {}: \t{}\n",i, time_sv_p.elapsed().as_nanos());

        let time_sv_v = ark_std::time::Instant::now();
        assert!(marlin_srs_verify(g1_generator, g2_generator, &srsj.srs, &pis, "verifier", c_size));
        print!("marlin_verification (verifier) {}: \t{}\n",i, time_sv_v.elapsed().as_nanos());
    }

    // batched timing
    let mut rand_r : Vec::<Vec::<E::Fr>> = Vec::new();
    let mut rand_t : Vec::<E::Fr> = Vec::new();
    rand_r.push(Vec::new());
    rand_r.push(Vec::new());
    rand_r.push(Vec::new());
    rand_r.push(Vec::new());
    for _ in 0..n_updates {
        rand_r[0].push(E::Fr::rand(rng));
        rand_r[1].push(E::Fr::rand(rng));
        rand_r[2].push(E::Fr::rand(rng));
        rand_r[3].push(E::Fr::rand(rng));
    }
    for _ in 0..2*c_size {
        rand_t.push(E::Fr::rand(rng));
    }
    for i in 0..n_iters {
        let time_svb_p = ark_std::time::Instant::now();
        assert!(marlin_srs_verify_batched(g1_generator, g2_generator, &srsj.srs, &pis, "prover",c_size, &rand_r, &rand_t));
        print!("marlin_verification (prover) {}: \t{}\n",i, time_svb_p.elapsed().as_nanos());

        let time_svb_v = ark_std::time::Instant::now();
        assert!(marlin_srs_verify_batched(g1_generator, g2_generator, &srsj.srs, &pis, "verifier", c_size, &rand_r, &rand_t));
        print!("marlin_verification (verifier) {}: \t{}\n",i, time_svb_v.elapsed().as_nanos());
    }
}

fn main() {
    let n_iters = 1;
    let n_updates = 5;
    let c_size = 50;
    
    print!("\ni = {}, N = {}\n", n_updates, c_size);
    marlin_timing::<Bls12_377>(n_iters, n_updates, c_size);
    
}
