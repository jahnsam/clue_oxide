## CluE
CluE (Cluester Evolution) is a quantum spin dynamics program to simulate central spin decoherence,  CluE implements Yang and Liu's cluster correlation expansion [1,2].

### Installation
[Rust](https://www.rust-lang.org/)

Before compiling CLuE, it is a good idea to check that everything is working properly by running the the tests. 
```
cargo test
```
To build Clue, run the following
```
cargo build --release
```
### Options: #[config]
| field        | values| description |
|--------------|:-----:|-----------|
clash\_distance\_pbc | float | minimum distance in Å for which particle in periodic boundary condition copies are allowed to overlap particle in the primary cell|
|  cluster\_batch\_size | int | number of clusters to evaluate between saves |
| cluster\_method | r2cce, cce | method of evaluating clusters |
|  detected\_spin\_position | grid | coordinates in Å of the detected spin|
|  input\_structure\_file | string | path to the input PDB |
| load\_geometry | cube, sphere | shape of the spin system |
| magnetic\_field | float | magnetic field strength in T of the applied magnetic field|
|  max\_cluster\_size | int | maximum cluster size to be evaluated |
|  neighbor\_cutoff\_3\_spin\_hahn\_mod\_depth | float | minimum modulation depth for the three spin hahn echo for an edge to be placed between two spins|
|neighbor\_cutoff\_3\_spin\_hahn\_taylor\_4 | float |  minimum magnitude of the fourth order coefficient, in the three spin Hahn echo Taylor series, (expanded in time),  for an edge to be placed between two spins |
| neighbor\_cutoff\_delta\_hyperfine | float | minimum magnitude of differnece of the _zz_-hyperfine in Hz, for an edge to be placed between two spins|
| neighbor\_cutoff\_dipole\_dipole | float |  minimum absolute __zz__-dipole-dipole coupling in Hz for an edge to be placed between two spins|
| neighbor\_cutoff\_dipole\_perpendicular | float | minimum absolute perpendicular-dipole-dipole coupling in Hz for an edge to be placed between two spins|
| number\_timepoints  | [int] |  number of timepoints used for each time increment| 
| pulse\_sequence | cp-_n_, hahn | pulse sequence to simulate with _n_ π-pulses|
| radius | float | radius in Å from the average position of the detected spin to use|
|  temperature | float | temperature in K used to calculate each clusters density matrix|
|  time\_increments  | [float] |  list of time increments in s used for propagation| 
|  write\_auxiliary\_signals | string |  file name to save cluster auxiliary signals under|
|  write\_bath | string | file name to save bath spins under|
|  write\_clusters | string | file name to save clusters under|
|  write\_info | string |  directory name to save general information under|
|  write\_exchange\_groups | string |  file name to save exchange groups (such as methyls) under|
|  write\_structure\_pdb | string |  PDB file name to save the spin system under|

## pyCluE
```
cd pyclue
```

```
maturin build
``` 
## References
<a id="1">[1]</a> 
Yang, W.; Liu, R. B. Decoherence of Coupled Electron Spins via Nuclear Spin Dynamics in Quantum Dots. Phys. Rev. B 2008, 77 (8), 085302. https://doi.org/10.1103/PhysRevB.77.085302.

<a id="2">[2]</a> 
Yang, W.; Liu, R.-B. Quantum Many-Body Theory of Qubit Decoherence in a Finite-Size Spin Bath. II. Ensemble Dynamics. Phys. Rev. B 2009, 79 (11), 115320. https://doi.org/10.1103/PhysRevB.79.115320.

<a id="2">[3]</a> 
Witzel, W. M.; Das Sarma, S. Quantum Theory for Electron Spin Decoherence Induced by Nuclear Spin Dynamics in Semiconductor Quantum Computer Architectures: Spectral Diffusion of Localized Electron Spins in the Nuclear Solid-State Environment. Phys. Rev. B 2006, 74 (3), 035322. https://doi.org/10.1103/PhysRevB.74.035322.

