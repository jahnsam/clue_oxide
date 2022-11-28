#[define spins]
// format: index particle
1 1H;
2 14N;
3 1H;


#[tensors]
// format: indices tesor
//
// vec^T * spin_op_vec: 
// index, [array of 3 elements] 
0, [0,0,ez]; // index 0 is reserved for the detected spin

1, [0,0,nz1];
2, [0,0,nz2];

// Duplicate entries add.
2, [Bx,0,0];

// spin_op_vec^T * matrix * spin_op_vec: 
// idx1 idx2 [array of 1, 3, or 9 elements] 
//
// Use 9 indices to specify the full tensor.
//

0, 1, [A11, A12, A13
       A21, A22, A23,
       A31, A32, A33];

// Use 3 indices to specify a diagonal matrix.
0, 2, [T11, T22, T33];

// Use 1 index to specify a matrix proportional to the identity.
0, 2, [fc];


// section header for the clusters
#[clusters]
1, 2;
1, 3;
2, 3;
1, 2, 3;


