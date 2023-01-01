use crate::space_3d::Vector3D;

use lebedev_laikov;

pub struct IntegrationGrid{
  dim: usize,
  points: Vec::<f64>,
  weights: Vec::<f64>,
}

impl IntegrationGrid{

  pub fn lebedev(n: usize) -> Self{
    let (x, y, z, weights) = lebedev_laikov::ld_vecs(n);

    let dim = 3;
    let n_elements = n*dim;

    let mut points = Vec::<f64>::with_capacity(n_elements);

    for ((xi,yi),zi) in x.iter().zip(y.iter()).zip(z.iter()){
      points.push(*xi);
      points.push(*yi);
      points.push(*zi);

    }

    IntegrationGrid{dim,points,weights}

  }

  //----------------------------------------------------------------------------
  pub fn len(&self) -> usize { self.weights.len() }
  //----------------------------------------------------------------------------
  pub fn weight(&self, index: usize) -> f64{ 
    self.weights[index]
  }
  //----------------------------------------------------------------------------
  pub fn x(&self, index: usize) -> f64{ 
    let idx = self.dim*index;
    self.points[idx]
  }
  //----------------------------------------------------------------------------
  pub fn y(&self, index: usize) -> f64{ 
    let idx = self.dim*index;
    self.points[idx+1]
  }
  //----------------------------------------------------------------------------
  pub fn z(&self, index: usize) -> f64{ 
    let idx = self.dim*index;
    self.points[idx+2]
  }
  //----------------------------------------------------------------------------
  pub fn xyz(&self, index: usize) -> Vector3D{
    Vector3D::from( [ self.x(index), self.y(index), self.z(index) ] )
  }
  //----------------------------------------------------------------------------

  pub fn remove_3d_hemisphere(self) -> IntegrationGrid{

    let mut points = Vec::<f64>::with_capacity(self.points.len()/2);
    let mut weights = Vec::<f64>::with_capacity(self.len()/2);

    let thr = 1e-12;
    for ii in 0..self.len(){
      let x = self.x(ii);
      let y = self.y(ii);
      let z = self.z(ii); 
      
      if (z < 0.0) || (z < thr && y < 0.0) || (z < thr && y < thr && x < 0.0){ 
        continue; 
      }

      points.push(x);
      points.push(y);
      points.push(z);
      weights.push(self.weight(ii));
    }

    let norm: f64 = weights.iter().sum();

    IntegrationGrid{
    dim: self.dim,
    points,
    weights: weights.iter().map(|x| x/norm).collect()
    }
  
  }
  //----------------------------------------------------------------------------
 
}



#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_lebedev(){
    let n = 6;
    let grid = IntegrationGrid::lebedev(6);
    let norm: f64 = grid.weights.iter().sum();
    assert!( (norm- 1.0).abs() <1e12);

    assert_eq!(grid.points, vec![
        1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, -1.0,
    ])
  } 
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_3d_hemisphere(){
    let n = 6;
    let grid = IntegrationGrid::lebedev(6).remove_3d_hemisphere();
    let norm: f64 = grid.weights.iter().sum();
    assert!( (norm- 1.0).abs() <1e12);
    assert_eq!(grid.points, vec![
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    ])

  } 
}






