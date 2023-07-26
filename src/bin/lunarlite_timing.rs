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
use uu_snarks::lunarlite_srs::*;
use ark_bls12_377::Bls12_377;

fn lunarlite_timing<E: PairingEngine>(n_iters:usize, n_updates: usize, c_size: usize)
{
    let rng = &mut test_rng();
    
    let mut pis: Vec<LunPi<E>> = Vec::new();
    let time_sg = ark_std::time::Instant::now();
    let mut srsj = lunarlite_srs_gen(rng, c_size);
    print!("lunarlite_gen: \t{}\n", time_sg.elapsed().as_nanos());
    pis.push(srsj.pi);

    for _ in 1..n_updates {
        let time_su = ark_std::time::Instant::now();
        srsj = lunarlite_srs_update(rng,&srsj.srs,&pis,c_size);
        print!("lunarlite_update: \t{}\n", time_su.elapsed().as_nanos());
        pis.push(srsj.pi);
    }
    
    // regular timing
    for i in 0..n_iters {
        let time_sv_p = ark_std::time::Instant::now();
        assert!(lunarlite_srs_verify(&srsj.srs, &pis, "prover", c_size));
        print!("lunarlite_verification (prover) {}: \t{}\n",i, time_sv_p.elapsed().as_nanos());

        let time_sv_v = ark_std::time::Instant::now();
        assert!(lunarlite_srs_verify(&srsj.srs, &pis, "verifier", c_size));
        print!("lunarlite_verification (verifier) {}: \t{}\n",i, time_sv_v.elapsed().as_nanos());
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
    assert!(lunarlite_srs_verify_batched(&srsj.srs, &pis, "prover",c_size, &rand_st, &rand_h));
    print!("lunarlite_verification (prover): \t{}\n", time_svb_p.elapsed().as_nanos());

    let time_svb_v = ark_std::time::Instant::now();
    assert!(lunarlite_srs_verify_batched(&srsj.srs, &pis, "verifier", c_size, &rand_st, &rand_h));
    print!("lunarlite_verification (verifier): \t{}\n", time_svb_v.elapsed().as_nanos());
}

fn main() {
    let n_iters = 1;
    let n_updates = 5;
    let c_size = 50;
    
    print!("\ni = {}, N = {}", n_updates, c_size);
    lunarlite_timing::<Bls12_377>(n_iters, n_updates, c_size);
    
}
