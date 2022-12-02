//! CluE (Cluster Evolution) is a spin dynamics simulation program for
//! electron spin decoherence.
//!
//! CluE implements Yang an Liu's cluster correlation expansion for an input
//! structure, usually a PDB file.  Simulation details are specified in a config
//! file.
pub mod config;
pub mod exchange_groups;
pub mod clue_errors;
pub mod pdb;
pub mod particle_config;
pub mod particle_specifier;
pub mod physical_constants;
pub mod primary_structure;
pub mod vector3;

/// This is the function called when using CluE through the command line.
pub fn clue() {
println!("\n\n{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}",
"##############################################################################\n",
"#                                                                            #\n",
"#                      ______    _             ________                      #\n",
"#                     /  ____|  | |           | _______|                     #\n",
"#                    |  |       | |   _   _   | |_____                       #\n",
"#                    |  |       | |  | | | |  |  _____|                      #\n",
"#                    |  |____   | |  | |_| |  | |______                      #\n",
"#                     \\______|  |_|  \\_____/  |________|                     #\n",
"#                                                                            #\n",
"#                 ___               ______               ___                 #\n",   
"#                    \\             /      \\             /                    #\n",
"#                     \\           /        \\           /                     #\n",
"#                      \\    O    /          \\    O    /                      #\n",
"#                       \\       /            \\       /                       #\n",
"#                        \\     /              \\     /                        #\n",
"#                ()       \\   /       ()       \\   /       ()                #\n",
"#                          \\ /                  \\ /                          #\n",
"#                           V                    V                           #\n",
"#                                                                            #\n",
"##############################################################################\n");
    println!("CLuE Copyright (C) 2022 Samuel M. Jahn");
    println!("This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.");
    println!("This is free software, and you are welcome to redistribute it");
    println!("under certain conditions; type `show c' for details.\n");


}
