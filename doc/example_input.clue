
// comment

/*
   another comment
*/

// section header for general configurations
// optional when it is the first/only section.
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


#[filter(label=water_h)]
  element in [H];
  residue in [WAT,SOL]
  
#[filter(label=glycerol_oh)]
  element in [H];
  residue in [GLYOL]  
  bonded_elements in [O]

#[filter(label=glycerol_ch)]
  element in [H];
  residue in [GLYOL]  
  bonded_elements in [C]

#[filter(label=sol_o)]
  element in [O];
  residue in [WAT,SOL,GLYOL]
  

#[filter(label=methyl_h)]
  element in [H];
  exchange_group in [methyl]


#[filter(label=TEMPO_N)]
  element in [N];
  
#[filter(label=TEMPO_axis)]
  indices in [17,18];


#[property(label= water_h)]
// for diff(a,b), both a and b are secondary filters 
// and must resovle to a single partilce.
// 
// The keyword "particle" refers to same does the property being set.
//   
// The rhs of axis properties should resolve to a [f64; 3].
//  
hf_z_axis = diff(particle, connected(sol_o) );  
hf_z_axis = diff(particle, connected(sol_o) );  
hf_xy_axis = diff(particle, same_molecule(water_h) );  

#[property(label = methyl_h)]
tunnel_splitting = 80e3; // Hz

#[property(label = TEMPO_N)]
hf_xy_axis = vector(TEMPO_axis); // TODO: how should this work?
hf_xy_axis = [1,0, 0.0 ,0.0]; 

#[property(label = glycerol_oh)]
nq_z_axis= diff(particle, connected(sol_o));  

#[property(label = glycerol_oh)]
cosubstitute = same_residue_sequence(particle);

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
F[0] =  [0,0,ez]; // index 0 is reserved for the detected spin

F[1] = [0,0,nz1];
F[2] = [0,0,nz2];

// entries can be added to
F[2] += [Bx,0,0];

// spin_op_vec^T * matrix * spin_op_vec: 
// idx1 idx2 [array of 1, 3, or 9 elements] 
//
// Use 9 numbers to specify the full tensor.
//

T[0,1] = [A11, A12, A13
     A21, A22, A23,
     A31, A32, A33];

// Use 6 number to specify a symmetric tensor.
T[0,3] = [A11, A12, A13
               A22, A23,
                    A33];

// Use 3 indices to specify a diagonal matrix.
T[0,2] = [T11, T22, T33];

// Use 1 index to specify a matrix proportional to the identity.
T[0,2] = [fc];


// section header for the clusters
#[clusters]
[1, 2];
[1, 3];
[2, 3];
[1, 2, 3];


