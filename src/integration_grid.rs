use crate::space_3d::Vector3D;

use lebedev_laikov;

/// `IntegrationGrid` is a struct for use in 
/// [quadrature](https://en.wikipedia.org/wiki/Quadrature_(mathematics)).
#[derive(Debug, Clone,PartialEq)]
pub struct IntegrationGrid{
  dim: usize,
  points: Vec::<f64>,
  weights: Vec::<f64>,
}

impl IntegrationGrid{

  /// This function generates a 
  /// [Lebedev](https://en.wikipedia.org/wiki/Lebedev_quadrature)
  /// quadrature grid as an `IntegrationGrid`.
  /// It uses the
  /// [lebedev_laikov crate](https://github.com/Rufflewind/lebedev_laikov)
  /// by [Rufflewind](https://github.com/Rufflewind).
  /// # Example
  /// ```
  ///# use clue_oxide::integration_grid::IntegrationGrid; 
  /// let n = 26;
  /// let grid = IntegrationGrid::lebedev(n);
  /// ```
  ///
  /// # Panic
  /// `lebedev()` will paninc if
  /// `n` ∉  \{6,   14,   26,   38,   50,   74,   86,  110,  146,  170,
  ///        194,  230,  266,  302,  350,  434,  590,  770,  974, 1202,
  ///       1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,
  ///       5294, 5810\}.
  ///
  pub fn lebedev(n: usize) -> Self{

    const LEBEDEVGRIDSIZES: [usize;32] = [
               6,   14,   26,   38,   50,   74,   86,  110,  146,  170, 
             194,  230,  266,  302,  350,  434,  590,  770,  974, 1202, 
            1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 
            5294, 5810];

    assert!(LEBEDEVGRIDSIZES.contains(&n));

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
  pub fn new(dim: usize) -> Self {
    let points = Vec::<f64>::new();
    let weights = Vec::<f64>::new();

    IntegrationGrid{dim,points,weights}
  }
  //----------------------------------------------------------------------------
  pub fn push(&mut self, points: Vec::<f64>,weight: f64){
    assert_eq!(points.len(), self.dim);

    for x in points.into_iter(){
      self.points.push(x);
    }

    self.weights.push(weight);
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
    assert!(self.dim >= 2);
    let idx = self.dim*index;
    self.points[idx+1]
  }
  //----------------------------------------------------------------------------
  pub fn z(&self, index: usize) -> f64{ 
    assert!(self.dim >= 3);
    let idx = self.dim*index;
    self.points[idx+2]
  }
  //----------------------------------------------------------------------------
  pub fn xyz(&self, index: usize) -> Vector3D{
    assert_eq!(self.dim,3);
    Vector3D::from( [ self.x(index), self.y(index), self.z(index) ] )
  }
  //----------------------------------------------------------------------------
  pub fn mean(&self) -> Vec::<f64>{
    let mut r: Vec::<f64> = (0..self.dim).map(|_| 0.0).collect();

    for (ii,x) in self.points.iter().enumerate(){
      let idx = ii%self.dim;
      r[idx] += x;
    }

    let normalizing_factor = (self.dim as f64)/(self.points.len() as f64);
    r = r.iter().map(|x| normalizing_factor*x).collect();

    r
  }
  //----------------------------------------------------------------------------
  pub fn translate(&mut self, vector: &Vec::<f64>){
    assert_eq!(vector.len(), self.dim);

    for (ii,x) in self.points.iter_mut().enumerate(){
      let idx = ii%self.dim;
      *x += vector[idx];
    }

  }
  //----------------------------------------------------------------------------

  pub fn remove_3d_hemisphere(self) -> IntegrationGrid{
    assert_eq!(self.dim,3);

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
    ]);

    assert_eq!(grid.mean(), vec![0.0,0.0,0.0]);
  } 
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_3d_hemisphere(){
    let grid = IntegrationGrid::lebedev(6).remove_3d_hemisphere();
    let norm: f64 = grid.weights.iter().sum();
    assert!( (norm- 1.0).abs() <1e12);
    assert_eq!(grid.points, vec![
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    ]);

    assert_eq!(grid.mean(), vec![1.0/3.0,1.0/3.0,1.0/3.0]);

  } 
  //----------------------------------------------------------------------------
  #[test]
  fn test_translate(){
    let mut grid = IntegrationGrid::new(3);
    let w = 1.0/6.0;
    grid.push(vec![1.0, 2.0, 2.0],w);
    grid.push(vec![3.0, 2.0, 2.0],w);
    grid.push(vec![2.0, 1.0, 2.0],w);
    grid.push(vec![2.0, 3.0, 2.0],w);
    grid.push(vec![2.0, 2.0, 1.0],w);
    grid.push(vec![2.0, 2.0, 3.0],w);
    assert_eq!(grid.mean(), vec![2.0, 2.0, 2.0]);
    let r = vec![-2.0,-2.0,-2.0];
    grid.translate(&r);
    assert_eq!(grid.mean(), vec![0.0, 0.0, 0.0]);
  }
}






