
// comment

/*
   another comment
*/

// section header for general configurations
#[config]
magnetic_field = 1.2; // T

pulse_sequence = hahn;



// section header for pdb info
#[pdb, path=path/to/pdb]

// If a path is put into a header, the section's contents can be placed
// in another file.  
// For example, a pdb can be used by specifiying the path to it
// or the pdb can be copy-pasted here

// section header for convergence search
//#[converge]
//system_radius
//neighbor_
//orientations 

// section header for convergence search
//#[cutoff(neighbor)]
//parameter = delta_hyperfine;
//#[neighbor_cutoff(delta_hyperfine)]
#[explore(delta_hyperfine), type = neighbor_cutoff]
value = 3e6; // Hz
range = [0,3e6]; // Hz
update = divide_by(10^0.2);
convergence_metric = rmsd;
threshold = 1e-3;
filter = water_h; // How does this work for neighbor cutoffs?
//filter_style = both;
//filter_style = neither;
filter_style = either; 


// section header for detected spin
#[detect]

// Specify the detected spin by a 3-vector.
//position = [x,y,z]
// or by pdb indices.
position = average_position_of_pdb_indices(28,29); 
spin = 1/2;
transition = [-1/2, +1/2];


// Define particle specific properties
// Each particle is covered covered exactly one header.
// In the case that a particle could fall under multiple headers,
// the ambiguity is resolved by assining it to the last header defined.
//
// CluE pre-defines a catch-all category that all particles default to
// unless a user defined category applies.
//
// Examples:
//
#[particle]
// Specifiers are preceded with the keyword "with".
// Specifiers can require the presence of some properties.
with element = N;
// Specifiers can also require the absence of some properties.
with residue != TEMPO;

// Specifiers can require that one ove several properties apply.
with pdb_ids in [1,2,3];

// Property setters are preceded with the keyword "set".
set isotope = 15N

set hyperfine_tensor = [Ax,Ay,Az];
set hyperfine_x_axis = vector_between_pdb_indices(1,19);
set hyperfine_y_axis = vector_between_pdb_indices(28,29);

set electric_quadrupole_tensor = [q1,q2,q3];
set hyperfine_x_axis = pdb_indices(1,19);

// This should err since the nitrogen in TEMPO is bonded to two carbons,/ 
// set electric_quadrupole_x_axis = connected(carbon, carbon);

// The nitrogen is only bonded to a single oxygen, so this should work.
set electric_quadrupole_y_axis = connected(self, oxygen);

#[particle]
with exchange_group in [methyl];
with connected_element in [C];
with residue in [TEMPO]
set tunnel_splitting = 80e3; // Hz

#[particle]
with exchange_group in [methyl];
with connected_element in [C];
with residue not in [TEMPO]
set tunnel_splitting = 100e3; // Hz

#[particle]
with residue in [WAT];
with element in [H];
set electric_quadrupole_tensor = [q1,q2,q3];
set electric_quadrupole_z_axis = connected(self, oxygen);

// This works because although water has two hydrogens,
// one of the hydrogens is the one for which the tensor applies 
// (specified by "self"), leaving only one other hydrogen in the molecule.
set electric_quadrupole_x_axis = same_molecule(self,hydrogen);

#[particle]
with residue in [GLY];
with element in [H];
set electric_quadrupole_tensor = [q1,q2,q3];
set electric_quadrupole_x_axis = random;
set electric_quadrupole_z_axis = connected(self, oxygen);


#[define spins]
// format: index particle
1 1H;
2 14N;
3 1H;


#[tensors]
// format: indices tesor
//
// vec^T * spin_op_vec: 
// index [array of 3 elements] 
0 [0,0,ez]; // index 0 is reserved for the detected spin

1 [0,0,nz1];
2 [0,0,nz2];

// duplicate entries add
2 [Bx,0,0]

// spin_op_vec^T * matrix * spin_op_vec: 
// idx1 idx2 [array of 1, 3, or 9 elements] 
//
// Use 9 indices to specify the full tensor.
//

0 1 [A11, A12, A13
     A21, A22, A23,
     A31, A32, A33];

// Use 3 indices to specify a diagonal matrix.
0 2 [T11, T22, T33];

// Use 1 index to specify a matrix proportional to the identity.
0 2 [fc];


// section header for the clusters
#[clusters]
1 2
1 3
2 3
1 2 3


