use crate::space_3d::Vector3D;

use rand::distributions::Uniform;
use rand_distr::Distribution;
use rand_chacha::ChaCha20Rng;


pub fn get_kmeans_partition(
    rng: &mut ChaCha20Rng,
    k: usize,
    data: &[Vector3D],
    ) -> Vec::<usize>
{

  let mut klusters = initialize_klusters(rng, k, &data);

  assign_data(&mut klusters, &data);

  move_kluster_centers(&mut klusters, &data);

  let mut klusters0: Vec::<Kluster>;

  loop{
    klusters0 = klusters.clone();

    assign_data(&mut klusters, &data);

    move_kluster_centers(&mut klusters, &data);


    if are_klusters_identical(&klusters,&klusters0){
      break;
    }
    
  }

  let mut element_to_block = (0..data.len()).map(|_| 0)
    .collect::<Vec::<usize>>();

  for (hh,kluster) in klusters.iter().enumerate(){
    for idx in kluster.elements.iter(){
      element_to_block[*idx] = hh;
    }
  }

  element_to_block
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// K-Means Clustering Cluster
// Since `Cluster` is already an important data structure in CluE Oxide,
// the k-means cluster is called `Kluster`.
#[derive(Debug,Clone)]
struct Kluster{
  center: Vector3D,
  elements: Vec::<usize>,
}

impl Kluster{
  //----------------------------------------------------------------------------
  fn random(
    rng: &mut ChaCha20Rng,
    min: &Vector3D,
    max: &Vector3D,
  ) -> Self
  {
    let mut center = Vector3D::zeros();

    let range = Uniform::new(0.0f64, 1.0);

    for (ii,r) in center.elements.iter_mut().enumerate(){
      let delta_x = max.elements[ii] -min.elements[ii];
      let random_number = range.sample(rng);
      *r = min.elements[ii] + delta_x * random_number; 
    }
    Self{
      center,
      elements: Vec::<usize>::new(),
    }
  }
  //----------------------------------------------------------------------------
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//------------------------------------------------------------------------------
fn are_klusters_identical(klusters: &[Kluster],klusters0: &[Kluster]) -> bool
{

  if klusters.len() != klusters0.len(){
    return false;
  }

  for (hh,kluster) in klusters.iter().enumerate(){
    if kluster.elements != klusters0[hh].elements{
      return false;
    }
  }

  true
}
//------------------------------------------------------------------------------
fn move_kluster_centers(klusters: &mut Vec::<Kluster>, data: &[Vector3D])
{
  for kluster in klusters.iter_mut(){
    if kluster.elements.is_empty(){
      continue;
    }
    let mut center = Vector3D::zeros();
    for idx in kluster.elements.iter(){
      center = &center + &data[*idx];
    }
    center = center.scale(1.0/(kluster.elements.len() as f64));

    kluster.center = center;
  }
}
//------------------------------------------------------------------------------
fn assign_data(klusters: &mut Vec::<Kluster>, data: &[Vector3D]){

  for kluster in klusters.iter_mut(){
    kluster.elements = Vec::<usize>::new();
  }

  for (ii,x) in data.iter().enumerate(){
    let mut min_index = 0;
    let mut min_distance_2 = f64::INFINITY;

    for (hh,kluster) in klusters.iter().enumerate(){
      let distance_2 = (x - &kluster.center).norm_squared();

      if distance_2 < min_distance_2{
        min_distance_2 = distance_2;
        min_index = hh;
      }
    }
    klusters[min_index].elements.push(ii);
  }
}
//------------------------------------------------------------------------------
fn initialize_klusters(
    rng: &mut ChaCha20Rng,
    k: usize,
    data: &[Vector3D],
    ) -> Vec::<Kluster>
{
  let ones = Vector3D::from([1.0,1.0,1.0]);
  let mut max = ones.scale(-f64::INFINITY);
  let mut min = ones.scale(f64::INFINITY);
  
  for r in data.iter(){
    for (ii,x) in r.elements.iter().enumerate(){
      max.elements[ii] = f64::max(*x,max.elements[ii]);
      min.elements[ii] = f64::min(*x,min.elements[ii]);
    }
  }

  let mut klusters = Vec::<Kluster>::with_capacity(k);
  for _ii in 0..k{
    klusters.push(Kluster::random(rng,&min,&max));
  }
   
  klusters
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#[cfg(test)]
mod tests{
  use super::*;
  use rand::SeedableRng;

  //----------------------------------------------------------------------------
  // TODO: Improve test so that it passes for all seed values.
  // Since k-means can settle into one of multiple minima, only some 
  // seeds give the expected result. 
  #[test]
  fn test_get_kmeans_partition(){
    let k = 4;
    let mut rng = ChaCha20Rng::seed_from_u64(2);
    let (data,_expected_klusters) = get_test_kluster_data();

    let elements_to_block = get_kmeans_partition(&mut rng,k,&data);
    assert_eq!(elements_to_block.len(),12); 

    assert_eq!(elements_to_block[0],elements_to_block[1]);
    assert_eq!(elements_to_block[0],elements_to_block[2]);

    assert_eq!(elements_to_block[3],elements_to_block[4]);
    assert_eq!(elements_to_block[3],elements_to_block[5]);

    assert_eq!(elements_to_block[6],elements_to_block[7]);
    assert_eq!(elements_to_block[6],elements_to_block[8]);

    assert_eq!(elements_to_block[9],elements_to_block[10]);
    assert_eq!(elements_to_block[9],elements_to_block[11]);

    assert_ne!(elements_to_block[0],elements_to_block[3]);
    assert_ne!(elements_to_block[0],elements_to_block[6]);
    assert_ne!(elements_to_block[0],elements_to_block[9]);

    assert_ne!(elements_to_block[3],elements_to_block[6]);
    assert_ne!(elements_to_block[3],elements_to_block[9]);

    assert_ne!(elements_to_block[6],elements_to_block[9]);
  
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_move_kluster_centers(){

    let (data,expected_klusters) = get_test_kluster_data();

    let mut klusters = vec![
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![0,1,2],
      },
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![3,4,5],
      },
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![6,7,8],
      },
      Kluster{
        center: Vector3D::zeros(),
        elements: vec![9,10,11],
      },
    ];

    move_kluster_centers(&mut klusters, &data);

    for (ii,kluster) in klusters.iter().enumerate(){
      let err = (&kluster.center - &expected_klusters[ii].center).norm();
      assert!( err < 1e-12 );
    }

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_assign_data_are_klusters_identical(){

    let (data,expected_klusters) = get_test_kluster_data();

    let mut klusters = vec![
      Kluster{
        center: (&data[0] + &(&data[1]+&data[2])).scale(1.0/3.0),
        elements: vec![0,1,2,3],
      },
      Kluster{
        center: (&data[3] + &(&data[4]+&data[5])).scale(1.0/3.0),
        elements: vec![4,5],
      },
      Kluster{
        center: (&data[6] + &(&data[7]+&data[8])).scale(1.0/3.0),
        elements: vec![9,10,11],
      },
      Kluster{
        center: (&data[9] + &(&data[10]+&data[11])).scale(1.0/3.0),
        elements: vec![6,7,8],
      },
    ];
    
    assert!(!are_klusters_identical(&klusters, &expected_klusters));

    assign_data(&mut klusters,&data);
    for (ii,kluster) in klusters.iter().enumerate(){
      assert_eq!(kluster.elements, expected_klusters[ii].elements)
    }

    assert!(are_klusters_identical(&klusters, &expected_klusters));

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_initialize_klusters(){
    let k = 4;
    let mut rng = ChaCha20Rng::from_entropy();
    let (data,_) = get_test_kluster_data();
    let klusters = initialize_klusters(&mut rng, k, &data);

    for kluster in klusters.iter(){
      let center = &kluster.center;
      assert!(center.x() > 34.530);
      assert!(center.x() < 38.970);
      assert!(center.y() > 34.140);
      assert!(center.y() < 40.260);
      assert!(center.z() > 35.140);
      assert!(center.z() < 39.500);
    }
  }
  //----------------------------------------------------------------------------
  fn get_test_kluster_data() -> (Vec::<Vector3D>,Vec::<Kluster>){  
    let data = vec![
      Vector3D::from([38.970,  35.510,  38.800]),
      Vector3D::from([37.380,  36.130,  39.500]),
      Vector3D::from([38.800,  37.180,  39.010]),

      Vector3D::from([37.110,  34.140,  37.860]),
      Vector3D::from([38.560,  34.210,  37.060]),
      Vector3D::from([37.060,  34.520,  36.030]),

      Vector3D::from([37.570,  38.870,  38.780]),
      Vector3D::from([35.840,  38.870,  38.940]),
      Vector3D::from([36.850,  40.260,  38.120]),

      Vector3D::from([35.260,  38.370,  35.140]),
      Vector3D::from([35.270,  39.920,  36.060]),
      Vector3D::from([34.530,  38.470,  36.820]),
    ];

    let klusters = vec![
      Kluster{
        center: (&data[0] + &(&data[1]+&data[2])).scale(1.0/3.0),
        elements: vec![0,1,2],
      },
      Kluster{
        center: (&data[3] + &(&data[4]+&data[5])).scale(1.0/3.0),
        elements: vec![3,4,5],
      },
      Kluster{
        center: (&data[6] + &(&data[7]+&data[8])).scale(1.0/3.0),
        elements: vec![6,7,8],
      },
      Kluster{
        center: (&data[9] + &(&data[10]+&data[11])).scale(1.0/3.0),
        elements: vec![9,10,11],
      },
    ];

    (data,klusters)
  }
  //----------------------------------------------------------------------------
}
