pub extern crate matrix_oxide;
pub use matrix_oxide as mox;
pub mod quantum;
pub mod tensors;
pub mod phys;
pub mod io;
pub mod Element;
pub mod Particle;
pub mod ParticleSpecifier;

//                  isotopologue                 orientation             CCE
// io -> config, structure -> spin_system -> clusters, tensors -> quantum -> signal -> io
//                                | ^                | ^  
//                                V |                V |
//                                io                 io

