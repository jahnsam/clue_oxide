//use crate::clue_errors::*;

use crate::cluster::{
  partition::PartitioningMethod,
  unit_of_clustering::UnitOfClustering,
};
use crate::config::*;
use crate::config::parse_properties::set_symmetric_tensor_3d;
use crate::physical_constants::ANGSTROM;
use crate::space_3d::Vector3D;
use crate::integration_grid::IntegrationGrid;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Config{
  //----------------------------------------------------------------------------
  /// This function parses a line of tokens and modifies the `Config`. 
  pub fn parse_config_line(&mut self, expression: &TokenExpression) 
    -> Result<(),CluEError>
  {

    if expression.relationship != Some(Token::Equals){
      return Err(CluEError::ExpectedEquality(expression.line_number));
    }

  let already_set = ||{
    CluEError::OptionAlreadySet(
        expression.line_number, expression.lhs[0].to_string()) };

    match expression.lhs[0]{
      Token::ApplyPBC 
        => set_to_some_bool(&mut self.apply_pbc,expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClashDistancePBC => {
        set_to_some_f64(&mut self.clash_distance_pbc, expression)?;
        if let Some(r) = &mut self.clash_distance_pbc{
            *r *= ANGSTROM;
        }else{ return Err(CluEError::NoClashDistancePBC);}
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClusterBatchSize 
        => set_to_some_usize(&mut self.cluster_batch_size, expression)?, 
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClusterDensityMatrix => {
        if let Some(_value) = &self.density_matrix{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }

        match rhs[0]{
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ApproxThermal | Token::Thermal => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let new_expression = TokenExpression{
              lhs: expression.lhs.clone(),
              rhs: Some(args),
              relationship: Some(Token::Equals),
              line_number: expression.line_number
            };

            let mut temperature_opt: Option<f64> = None;
            set_to_some_f64(&mut temperature_opt,&new_expression)?;
            if let Some(temperature) = temperature_opt{
              match rhs[0]{
                Token::ApproxThermal => self.density_matrix 
                  = Some(DensityMatrixMethod::ApproxThermal(temperature)),
                Token::Thermal => self.density_matrix 
                  = Some(DensityMatrixMethod::Thermal(temperature)),
                _ => return Err(CluEError::InvalidToken(expression.line_number,
                  expression.lhs[0].to_string())),
              }
            }else{
              return Err(CluEError::NoDensityMatrixMethod);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::Identity => self.density_matrix 
            = Some(DensityMatrixMethod::Identity),
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClusterMethod 
        => set_to_some_cluster_method(&mut self.cluster_method,expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::ClustersFile 
        => set_to_some_string(&mut self.clusters_file, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::DetectedSpinGMatrix 
        | Token::DetectedSpinGX| Token::DetectedSpinGY | Token::DetectedSpinGZ
        => set_symmetric_tensor_3d(&mut self.detected_spin_g_matrix,
            &expression.lhs[0], expression, "detected_spin")?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::DetectedSpinPosition =>{
        if let Some(_value) = &self.detected_spin_position{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }

        match rhs[0]{
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::Centroid =>{
            let args = get_function_arguments(rhs,expression.line_number)?;
            if args.is_empty(){
              return Err(CluEError::InvalidArgument(expression.line_number,
                    "#[group]".to_string()));
            }

            let filter_labels = vec_tokens_to_vec_strings(args)?;

            self.detected_spin_position = Some(
                DetectedSpinCoordinates::CentroidOverGroup(filter_labels));

          }
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::CentroidOverSerials =>{
            let args = get_function_arguments(rhs,expression.line_number)?;

            let mut serials = Vec::<u32>::new();
            let value_token = to_i32_token(args, expression.line_number)?;
            match value_token {
              Token::Int(a) => serials.push(a as u32),
              Token::VectorI32(v) => serials = v.iter().map(|&a| a as u32)
                .collect(),  
              _  => return Err(CluEError::InvalidArgument(expression.line_number,
                    "list of positive integers".to_string())),
            }
            self.detected_spin_position = Some(
                DetectedSpinCoordinates::CentroidOverSerials(serials));

          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ReadCSV => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let grid = IntegrationGrid::read_from_csv(&(args[0].to_string()))?
              .scale(ANGSTROM);

            if grid.dim() != 3{
              return Err(CluEError::WrongProbabilityDistributionDim(
                    expression.line_number,3,grid.dim()));
            }
            self.detected_spin_position = Some(
                DetectedSpinCoordinates::ProbabilityDistribution(grid));

          }
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::SquareBracketOpen =>{
            let mut r_opt: Option<Vector3D> = None;
            set_to_some_vector3d(&mut r_opt, expression)?;
            if let Some(r_angstrom) = r_opt{
              let r = r_angstrom.scale(ANGSTROM); 
              self.detected_spin_position 
                = Some(DetectedSpinCoordinates::XYZ(r));
            }else{
              return Err(CluEError::NoCentralSpinCoor);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
        }
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::InputStructureFile 
        => set_to_some_string(&mut self.input_structure_file, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::LoadGeometry =>{
        if let Some(_value) = &self.load_geometry{
          return Err(already_set());
        }

        let Some(rhs) = &expression.rhs else { 
          return Err(CluEError::NoRHS(expression.line_number));
        };

        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }else if rhs.len() > 1{
          return Err(CluEError::TooManyRHSArguments(expression.line_number));
        }

        match rhs[0]{
          Token::Cube => self.load_geometry = Some(LoadGeometry::Cube),
          Token::Ball => self.load_geometry = Some(LoadGeometry::Ball),
          _ => return Err(CluEError::InvalidGeometry(expression.line_number,
                rhs[0].to_string())),
        }

      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MagneticField => {
  
        if let Some(_value) = &self.magnetic_field{
          return Err(already_set());
        }
        let mut mf_opt: Option<f64> = None;
        set_to_some_f64(&mut mf_opt, expression)?;

        if let Some(bz) = mf_opt{
          self.magnetic_field = Some( Vector3D::from([0.0,0.0,bz]) );
        }else{
          return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string()))
        }

      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MaxCellSize 
        => set_to_some_usize(&mut self.max_cell_size, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MaxClusterSize 
        => set_to_some_usize(&mut self.max_cluster_size, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MaxSpinOrder 
        => set_to_some_usize(&mut self.max_spin_order, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::MinCellSize 
        => set_to_some_usize(&mut self.min_cell_size, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffCoupling 
        => set_to_some_f64(&mut self.neighbor_cutoff_coupling,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDeltaHyperfine 
        => set_to_some_f64(&mut self.neighbor_cutoff_delta_hyperfine,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDipolePerpendicular
        => set_to_some_f64(&mut self.neighbor_cutoff_dipole_perpendicular,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoffDistance 
        => {
          set_to_some_f64(&mut self.neighbor_cutoff_distance,expression)?;
          if let Some(r) = &mut self.neighbor_cutoff_distance{
            *r *= ANGSTROM;
          }else{
            return Err(CluEError::NoNeighborCutoffDistance);
          }
        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoff3SpinHahnModDepth 
        => set_to_some_f64(&mut self.neighbor_cutoff_3_spin_hahn_mod_depth,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NeighborCutoff3SpinHahnTaylor4 
        => set_to_some_f64(&mut self.neighbor_cutoff_3_spin_hahn_taylor_4,
            expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NumberSystemInstances
        => set_to_some_usize(&mut self.number_system_instances, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::NumberTimepoints 
        => set_to_vec_usize(&mut self.number_timepoints, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::OrientationGrid => {
        if let Some(_value) = &self.orientation_grid{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }

        match rhs[0]{
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::Lebedev | Token::Random => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let new_expression = TokenExpression{
              lhs: expression.lhs.clone(),
              rhs: Some(args),
              relationship: Some(Token::Equals),
              line_number: expression.line_number
            };

            let mut n_grid_opt: Option<usize> = None;
            set_to_some_usize(&mut n_grid_opt,&new_expression)?;
            if let Some(n_grid) = n_grid_opt{
              match rhs[0]{
                Token::Lebedev => self.orientation_grid 
                  = Some(OrientationAveraging::Lebedev(n_grid)),
                Token::Random => self.orientation_grid 
                  = Some(OrientationAveraging::Random(n_grid)),
                _ => return Err(CluEError::InvalidToken(expression.line_number,
                  expression.lhs[0].to_string())),
              }
            }else{
              return Err(CluEError::NoOrientationGrid);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::ReadCSV => {
            let args = get_function_arguments(rhs,expression.line_number)?;

            if args.is_empty(){
              return Err(CluEError::TooFewRHSArguments(expression.line_number));
            }else if args.len() > 1{
              return Err(CluEError::TooManyRHSArguments(expression.line_number));
            }

            let grid = IntegrationGrid::read_from_csv(&(args[0].to_string()))?;

            if grid.dim() != 3{
              return Err(CluEError::WrongOrientationGridDim(
                    expression.line_number,3,grid.dim()));
            }
            self.orientation_grid = Some(
                OrientationAveraging::Grid(grid));

          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          Token::SquareBracketOpen =>{
            let mut r_opt: Option<Vector3D> = None;
            set_to_some_vector3d(&mut r_opt, expression)?;
            if let Some(r) = r_opt{
              let mut grid = IntegrationGrid::new(3);
              grid.push(vec![r.x(),r.y(),r.z()],1.0)?;
              self.orientation_grid
                = Some(OrientationAveraging::Grid(grid));
            }else{
              return Err(CluEError::NoOrientationGrid);
            }
          },
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
          _ => return Err(CluEError::InvalidToken(expression.line_number,
                rhs[0].to_string())),
          // -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::PartitioningMethod => {
        if let Some(_value) = &self.partitioning_method{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }
        self.partitioning_method 
            = Some(PartitioningMethod::from(rhs,expression.line_number)?);
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::PulseSequence => {
        if let Some(_value) = &self.pulse_sequence{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }
        match rhs[0]{
          Token::Hahn 
            => self.pulse_sequence = Some(PulseSequence::CarrPurcell(1)),
          Token::CarrPurcell => {
            if rhs.len() != 3 || rhs[1] != Token::Minus{
              return Err(CluEError::InvalidPulseSequence(
                    expression.line_number));
            }
            let rhs = token_stream::read_strings_as_integers(rhs.clone(), 
                expression.line_number)?;
            if let Token::Int(n_pi) = rhs[2]{
              self.pulse_sequence 
                = Some(PulseSequence::CarrPurcell(n_pi as usize));
            }else{
              return Err(CluEError::InvalidPulseSequence(
                    expression.line_number));
            }

          },
          _ => return Err(CluEError::InvalidPulseSequence(
                expression.line_number)),
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Radius => {
        set_to_some_f64(&mut self.radius,expression)?;
        if let Some(r) = &mut self.radius{
            *r *= ANGSTROM;
        }else{ return Err(CluEError::NoRadius);}
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::SaveDir 
        => set_to_some_string(&mut self.save_name, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::SystemName 
        => set_to_some_string(&mut self.system_name, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::RNGSeed => {
        if let Some(_value) = &self.rng_seed{
          return Err(already_set());
        }
        let mut rng_seed: Option<i32> = None;
        set_to_some_i32(&mut rng_seed,expression)?;

        if let Some(seed) = rng_seed{
          self.rng_seed = Some( seed as u64);
        }else{
          return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string()))
        }
      },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::Temperature
        => {
          //set_to_some_f64(&mut self.temperature, expression)?
          return Err(CluEError::DeprecatedKeywordReplaced(
                expression.line_number,
                "temperature = T;".to_string(),
                "cluster_density_matrix = thermal(T);".to_string(),
                ))
        },
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::TimeIncrements 
        => set_to_vec_f64(&mut self.time_increments, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::UnitOfClustering => {
        if let Some(_value) = &self.unit_of_clustering{
          return Err(already_set());
        }
        let Some(rhs) = &expression.rhs else{
          return Err(CluEError::NoRHS(expression.line_number));
        }; 
        if rhs.is_empty(){
          return Err(CluEError::NoRHS(expression.line_number));
        }
        self.unit_of_clustering 
            = Some(UnitOfClustering::from(rhs,expression.line_number)?);
      }
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteAuxiliarySignals
        => set_to_some_bool(&mut self.write_auxiliary_signals, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteBath
        => set_to_some_bool(&mut self.write_bath, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteClusters
        => set_to_some_bool(&mut self.write_clusters, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteExchangeGroups
        => set_to_some_bool(&mut self.write_exchange_groups, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteInfo
        => set_to_some_bool(&mut self.write_info, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteMethylPartitions
        => set_to_some_bool(&mut self.write_methyl_partitions, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteOrientationSignals
        => set_to_some_bool(&mut self.write_orientation_signals, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteSansSpinSignals 
        => set_to_some_bool(&mut self.write_sans_spin_signals, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteStructurePDB 
        => set_to_some_bool(&mut self.write_structure_pdb, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Token::WriteTensors 
        => set_to_some_bool(&mut self.write_tensors, expression)?,
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      _ => return Err(CluEError::InvalidToken(expression.line_number,
            expression.lhs[0].to_string())),
    }
    Ok(())
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use crate::structure::particle_filter::{
    SecondaryParticleFilter,VectorSpecifier};

  #[test]
  fn test_parse_config_line(){
    let expressions = get_tokens_from_line("\
        #[config]\n
        apply_pbc = false;
        clash_distance_pbc = 0.1;
        cluster_batch_size = 20000;
        cluster_density_matrix = thermal(20);
        cluster_method = cce;
        detected_spin_g_matrix = [2.0097, 2.0064, 2.0025];
        detected_spin_g_y = [-1.1500, -0.4700, 0.7100];
        detected_spin_g_x = diff(group(tempo_c1) , filter(tempo_c19) );
        detected_spin_position = centroid_over_serials([28,29]);
        input_clusters_file = \"clusters_file.txt\";
        input_structure_file = \"../../assets/TEMPO_wat_gly_70A.pdb\";
        load_geometry = cube;
        max_cell_size = 2;
        max_cluster_size = 4;
        max_spin_order = 8;
        min_cell_size = 1;
        number_timepoints = [40,60];
        neighbor_cutoff_delta_hyperfine = 1e04;
        neighbor_cutoff_coupling = 1e+3;
        neighbor_cutoff_dipole_perpendicular = 100;
        neighbor_cutoff_3_spin_hahn_mod_depth = 1e-10;
        neighbor_cutoff_3_spin_hahn_taylor_4 = 1e-9;
        orientation_grid = lebedev(170);
        partitioning_method = particles;
        pulse_sequence = cp-1;
        radius = 80;
        save_dir = \"save_directory\";
        time_increments = [1e-9, 5e-7];
        write_auxiliary_signals = true;
        write_bath = true;
        write_clusters = true;
        write_exchange_groups = true;
        write_info = true;
        write_methyl_partitions = true;
        write_orientation_signals = true;
        write_sans_spin_signals = true;
        write_structure_pdb = true;
        write_tensors = false;
        ").unwrap();

    let mut config = Config::new();
    config.parse_token_stream(expressions).unwrap();
    config.set_defaults().unwrap();

    assert_eq!(config.apply_pbc, Some(false));
    let clash_distance_pbc = config.clash_distance_pbc.unwrap();
    assert!( (clash_distance_pbc - 0.1e-10).abs()/(0.1e-10) < 1e-12 );
    assert_eq!(config.cluster_batch_size,Some(20000));
    assert_eq!(config.density_matrix, Some(DensityMatrixMethod::Thermal(20.0)));
    assert_eq!(config.cluster_method,Some(ClusterMethod::CCE));

    let g_matrix = match config.detected_spin_g_matrix.unwrap(){
      TensorSpecifier::Eig(eig_specifier) => eig_specifier.clone(),
      _ => panic!("Expected TensorSpecifier::Eig(eig_specifier)."),
    };

    assert_eq!(g_matrix.values, 
        Some([2.0097, 2.0064, 2.0025] ) );
    assert_eq!(g_matrix.z_axis, None);
    assert_eq!(g_matrix.y_axis, 
        Some(VectorSpecifier::Vector(
            Vector3D::from([-1.1500, -0.4700, 0.7100]) ) ) );
    assert_eq!(g_matrix.x_axis, 
        Some(VectorSpecifier::Diff(
            SecondaryParticleFilter::Filter, "tempo_c1".to_string(),
            SecondaryParticleFilter::Filter, "tempo_c19".to_string(),
            )));

    assert_eq!(config.detected_spin_identity, Some(Isotope::Electron));
    assert_eq!(config.detected_spin_multiplicity, Some(2));
    assert_eq!(config.detected_spin_position, 
        Some(DetectedSpinCoordinates::CentroidOverSerials(vec![28,29])));
    assert_eq!(config.input_structure_file, 
         Some("../../assets/TEMPO_wat_gly_70A.pdb".to_string()));
    assert_eq!(config.clusters_file, 
         Some("clusters_file.txt".to_string()));
    
    assert_eq!(config.load_geometry, Some(LoadGeometry::Cube));
    assert_eq!(config.max_cell_size, Some(2));
    assert_eq!(config.max_cluster_size, Some(4));
    assert_eq!(config.max_spin_order, Some(8));
    assert_eq!(config.min_cell_size, Some(1));
    assert_eq!(config.neighbor_cutoff_delta_hyperfine, Some(1e4));
    assert_eq!(config.neighbor_cutoff_coupling, Some(1e3));
    assert_eq!(config.neighbor_cutoff_dipole_perpendicular, Some(1e2));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_mod_depth, Some(1e-10));
    assert_eq!(config.neighbor_cutoff_3_spin_hahn_taylor_4, Some(1e-9));
    assert_eq!(config.number_timepoints, vec![40,60]);
    assert_eq!(config.partitioning_method, Some(PartitioningMethod::Particles));
    assert_eq!(config.pulse_sequence, Some(PulseSequence::CarrPurcell(1)));
    assert_eq!(config.radius, Some(80.0e-10));
    assert_eq!(config.save_name, Some(String::from("save_directory")));
    assert_eq!(config.time_increments, vec![1e-9,5e-7]);
    assert_eq!(config.write_auxiliary_signals, Some(true));
    assert_eq!(config.write_bath, Some(true));
    assert_eq!(config.write_clusters, Some(true));
    assert_eq!(config.write_exchange_groups, Some(true));
    assert_eq!(config.write_info, Some(true));
    assert_eq!(config.write_methyl_partitions, Some(true));
    assert_eq!(config.write_orientation_signals, Some(true));
    assert_eq!(config.write_structure_pdb, Some(true));
    assert_eq!(config.write_sans_spin_signals, Some(true));
    assert_eq!(config.write_tensors, Some(false));
  }
  //----------------------------------------------------------------------------
}

