# BMS23: Benchmarking the Setup of Updatable zk-SNARKs

This code contains an implementation of the setup of updatable zk-SNARKs: SONIC, Marlin, Plonk, Lunar, and Basilisk. It can be used for benchmarking and comparing the performance of the setup phase in different schemes. 
The code is based on the [Arkworks](https://github.com/arkworks-rs) library. The description of the algorithms and results of this implementation can be found in the Latincrypt '23 paper, and its full version is available on the IACR eprint archive. The research was done by Karim Baghery, Axel Mertens, and Mahdi Sedaghat (Cosic, KU Leuven).

The structure of this repository is as follows:
* basilisk_srs.rs: Rust code implementing the setup of Basilisk.
* lunarlite_srs.rs: Rust code implementing the setup of LunarLite.
* marlin_srs.rs: Rust code implementing the setup of Marlin.
* sonic_srs.rs: Rust code implementing the setup of Sonic.
* bin: directory containing scripts for benchmarking:
  - basilisk_timing.rs: Rust code for benchmarking Basilisk.
  - identifiable_security.rs: Rust code for benchmarking the algorithm for identifiable security (applied to Basilisk).
  - lunarlite_timing.rs: Rust code for benchmarking LunarLite.
  - marlin_timing.rs: Rust code for benchmarking Marlin.
  - sonic_timing.rs: Rust code for benchmarking Sonic.

## How to use
The only prerequisites are rust and cargo to be installed. 

To build the code, and get all necessary packages from [Arkworks](https://github.com/arkworks-rs):
```
cargo build
```

After this, you can run eg. basilisk_timing.rs:
```
cargo run --bin basilisk_timing
```


## Example Benchmark
Say you want to benchmark Basilisk for a setting with 100 updates. First, set the upper bound on the size of SRS. For Basilisk, this is simply the number of multiplication gates in the circuit. In this example, let's say this is 50 000.

In ```main()``` of src/bin/basilisk_timing.rs, enter those values:
```
let n_updates = 100;
let c_size = 50000;
```
If you'd like to use a different curve than BLS12-381, you can edit this here as well. Make sure the curve is imported correctly.


Then, in the function ```basilisk_timing(...)```, change the protocol to what you want to happen. Make sure the randomness is produced outside of the timing, as this can be done offline. 
If we want to measure the time a verifier would need for verification of this setting, ```basilisk_timing(...)``` would look like this:
```
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

    let mut srsj = basilisk_srs_gen(g1_generator, g2_generator, rng, c_size);
    pis.push(srsj.pi);

    for i in 1..n_updates {
        srsj = basilisk_srs_update_parallel(g1_generator, g2_generator, rng, &srsj.srs, c_size);
        pis.push(srsj.pi);

        assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "prover", g1_generator, g2_generator,c_size, &rand_st, &rand_h));
        assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "verifier", g1_generator, g2_generator, c_size, &rand_st, &rand_h));
    }

    // batched timing
    let mut rand_st : Vec::<Vec::<E::Fr>> = Vec::new();
    let mut rand_h : Vec::<E::Fr> = Vec::new();
    for _ in 0..c_size {
        rand_h.push(E::Fr::rand(rng));
    }
    for _ in 0..n_updates {
        rand_st.push(vec!(E::Fr::rand(rng), E::Fr::rand(rng)));
    }

    let time_svb_v = ark_std::time::Instant::now();
    assert!(basilisk_srs_verify_batched(&srsj.srs, &pis, "verifier", g1_generator, g2_generator, c_size, &rand_st, &rand_h));
    print!("basilisk_verification (verifier): \t{}\n", time_svb_v.elapsed().as_nanos());

}
```
This function can be changed if you want to benchmark other parts of the setup. 

Now, run the code as described above, and the timings will be printed in the terminal.

Benchmarking other SNARKs works very similarly.
Note that you can use the implementation for Basilisk also for Plonk, only the calculation of the circuit size is different.

## License

This library is licensed under either of the following licenses, at your discretion.

 * Apache License Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

## Acknowledgements
This work has been supported in part by the Defense Advanced Research Projects Agency (DARPA) under contract No. HR001120C0085, by the FWO under an Odysseus project GOH9718N, by the Research Council KU Leuven C1 on Security and Privacy for Cyber-Physical Systems and the Internet of Things with contract number C16/15/058, and by CyberSecurity Research Flanders with reference number VR20192203.
