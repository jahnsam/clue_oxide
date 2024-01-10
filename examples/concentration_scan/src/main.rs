use clue_oxide as clue;
use clue_oxide::{
  physical_constants::Isotope,
  config::particle_config::{IsotopeAbundance,IsotopeDistribution}
};

fn main() {

    // Read input from a string rather than a file.
    let config = match clue::config::Config::from("
        input_structure_file = \"TEMPO_wat_gly_70A.pdb\";
        radius = 25; // angstroms.

        detected_spin_position = centroid_over_serials([28,29]);
        detected_spin_g_matrix = [2.0097, 2.0064, 2.0025];
        detected_spin_g_y = diff(group(tempo_n) , group(tempo_o) );
        detected_spin_g_x = diff(group(tempo_c1) , group(tempo_c19) );
        number_timepoints = [101,140];

        time_increments = [50,500]*1e-9;
        cluster_method = cce;
        max_cluster_size = 2;
        magnetic_field = 1.2; // T.
        pulse_sequence = hahn;

        number_system_instances = 10;
        orientation_grid = random(1);

        neighbor_cutoff_distance = 4; // angstroms.

      #[group(tempo_n)]
        residues in [TEM];
        elements in [N];

      #[group(tempo_o)]
        residues in [TEM];
        elements in [O];

      #[group(tempo_c1)]
        serials in [1];

      #[group(tempo_c19)]
        serials in [19];

      #[spin_properties(tempo_n, 14N)]
        hyperfine_coupling = [20,20,100]*1e6;
        hyperfine_x = diff(bonded(tempo_c1),bonded(tempo_c19));
        hyperfine_y = diff(particle,bonded(tempo_o));

        electric_quadrupole_coupling = [-0.28, -1.47, 1.75]*1e6;
        electric_quadrupole_x = diff(bonded(tempo_c1),bonded(tempo_c19));
        electric_quadrupole_y = diff(particle,bonded(tempo_o));

      #[group(solvent_h)]
        residues not in [TEM];
        elements in [H];

      #[structure_properties(solvent_h)]
        isotope_abundances = {1H: 1, 2H: 0};

    "){
      Ok(cfg) => cfg,
      Err(err) => {
        eprintln!("CluE Error: {}.",err);
        return;
      },
    };

 
    let mut solvent_h_idx: Option<usize> = None;

    for (idx, p_cfg) in config.particles.iter().enumerate(){
      if p_cfg.label == "solvent_h"{
        solvent_h_idx = Some(idx);
        break;
      }
    }

    let solvent_h_idx = solvent_h_idx.expect("Could not find solvent_h_idx.");


    // Initialize and loop over proton fractions.
    let proton_fractions = vec![1.0, 0.75, 0.5, 0.25, 0.0];
    for &frac in proton_fractions.iter(){
      
      let mut frac_config = config.clone();

      // Set save directory.
      frac_config.save_name = Some(format!("CluE-proton_frac_{}",frac));

      let properties = frac_config.particles[solvent_h_idx]
        .properties.as_mut().expect("No properties for solvent_h");

      let isotope_abundances = vec![
        IsotopeAbundance{isotope: Isotope::Hydrogen1, abundance: frac},
        IsotopeAbundance{isotope: Isotope::Hydrogen2, abundance: 1.0 - frac},
      ];

      properties.isotopic_distribution = IsotopeDistribution{
          isotope_abundances,
          void_probability: None
      };

      match clue::run(frac_config){
        Ok(_) => (),
        Err(err) => eprintln!("CluE Error: {}.",err),
      }
    }
}
