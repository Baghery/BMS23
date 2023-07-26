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
use uu_snarks::*;
use ark_ff::PrimeField;
use ark_ec::ProjectiveCurve;
use ark_bls12_377::Bls12_377;
use ark_bn254::Bn254;
use ark_bls12_381::Bls12_381;


fn sonic_timing<E: PairingEngine>(n_iters:usize, n_updates: usize, c_size: usize)
{
    let rng = &mut test_rng();
    
    let mut pis: Vec<SonPi<E>> = Vec::new();
    let time_sg = ark_std::time::Instant::now();
    let mut srsj: SonicSrs<E> = sonic_srs_gen(rng, c_size);
    print!("sonic_gen: \t{}\n", time_sg.elapsed().as_nanos());
    pis.push(srsj.pi);


    for _ in 1..n_updates {
        let time_su = ark_std::time::Instant::now();
        srsj = sonic_srs_update(rng, &srsj.srs, &pis, c_size);
        print!("sonic_update: \t{}\n", time_su.elapsed().as_nanos());
        pis.push(srsj.pi);
    }

    // regular timing 
    let time_sv_p = ark_std::time::Instant::now();
    assert!(sonic_srs_verify(&srsj.srs,&pis,"prover",c_size));
    print!("sonic_verify (prover): \t\t{}\n", time_sv_p.elapsed().as_nanos());
    let time_sv_v = ark_std::time::Instant::now();
    assert!(sonic_srs_verify(&srsj.srs,&pis,"verifier",c_size));
    print!("sonic_verify (verifier): \t\t{}\n", time_sv_v.elapsed().as_nanos());
 
    // batched timing
    let mut t = vec![vec![E::Fr::rand(rng)]];
    for _ in 0..2 {t.push(vec![E::Fr::rand(rng)])}
    for _ in 1..2*c_size {
        t[0].push(E::Fr::rand(rng));
        t[1].push(E::Fr::rand(rng));
    }
    let mut r: Vec<Vec<E::Fr>> = Vec::new();
    for _ in 0..4 {r.push(vec![E::Fr::rand(rng)])}
    for _ in 1..n_updates {
        r[0].push(E::Fr::rand(rng));
        r[1].push(E::Fr::rand(rng));
        r[2].push(E::Fr::rand(rng));
        r[3].push(E::Fr::rand(rng));
    }

    let time_svb_p = ark_std::time::Instant::now();
    assert!(sonic_srs_verify_batched(&srsj.srs, &pis, "prover", c_size, &r, &t));
    let bvtp = time_svb_p.elapsed().as_nanos();
    print!("sonic_verify_batched (prover): \t\t{}\n", bvtp);
    let time_svb_v = ark_std::time::Instant::now();
    assert!(sonic_srs_verify_batched(&srsj.srs, &pis, "verifier", c_size, &r, &t));
    let bvtv = time_svb_v.elapsed().as_nanos();
    print!("sonic_verify_batched (verifier): \t\t{}\n", bvtv);
}

fn main() {
    let n_iters = 1;
    let n_updates = 3;
    let c_size = 30;
    
    print!("\ni = {}, N = {}\n", n_updates, c_size);
    sonic_timing::<Bls12_377>(n_iters, n_updates, c_size);

}
