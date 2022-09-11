//extern crate matrix_oxide;
use matrix_oxide as mox;
//use num_complex::Complex;

//pub const I: Complex<f64> = Complex::<f64>{re: 0.0, im: 1.0};
//pub const ONE: Complex<f64> = Complex::<f64>{re: 1.0, im: 0.0};
pub const I: mox::Z64 = mox::I;
pub const ONE: mox::Z64 = mox::Z64{re: 1.0, im: 0.0};

// 2021: https://physics.nist.gov/cgi-bin/cuu/Value?hbar 
pub const HBAR: f64 = 1.054571817e-34; // Js.

// 2021: https://introcs.cs.princeton.edu/java/data/pi-10million.txt 
pub const PI: f64 = 3.141592653589793;
pub const ERROR_THRESHOLD: f64 = 1e-15;
