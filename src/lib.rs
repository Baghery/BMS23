#![cfg_attr(not(feature = "std"), no_std)]
#![warn(
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    missing_docs
)]
#![allow(clippy::many_single_char_names, clippy::op_ref, missing_docs, unused_mut)]
#![forbid(unsafe_code)]


#[macro_use]
extern crate ark_std;

#[cfg(feature = "r1cs")]
#[macro_use]
extern crate derivative;

pub mod basilisk_srs;
pub mod sonic_srs;
pub mod lunarlite_srs;
pub mod marlin_srs; 

pub use self::{basilisk_srs::*, sonic_srs::*, lunarlite_srs::*, marlin_srs::*};

use ark_std::vec::Vec;

