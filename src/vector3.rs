// 2021: https://introcs.cs.princeton.edu/java/data/pi-10million.txt
pub const PI: f64 = 3.141592653589793;

#[derive(Debug, Clone)]
pub struct Vector3{ 
  pub x: f64,
  pub y: f64,
  pub z: f64
}

#[derive(Debug, Clone)]
pub struct SphereVec3{ 
  pub r: f64,
  pub theta: f64,
  pub phi: f64
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
impl Vector3{

  pub fn new() -> Vector3{
    Vector3{
      x: 0.0,
      y: 0.0,
      z: 0.0,
    }
  }    

  pub fn from(r: [f64; 3]) -> Vector3{
    Vector3{
      x: r[0],
      y: r[1],
      z: r[2],
    }
  }

  pub fn magnitude(&self) -> f64 {                                              
    (self.x*self.x + self.y*self.y + self.z*self.z).sqrt()
  }                                                                                   
                                                                                    
  pub fn scale(&self,scale: f64) -> Vector3 {                                   
    Vector3{                                                              
      x: scale*self.x,
      y: scale*self.y,
      z: scale*self.z,
    }                                                                       
  }                                                                                
                                                                                 
  pub fn normalize(&self) -> Vector3 {                                        
    let r = self.magnitude();                                                        
    self.scale(1.0/r)                                                            
  }                                                                                
                                                                                    
  pub fn dot(&self, other: &Vector3) -> f64 {                                     
    self.x*other.x + self.y*other.y + self.z*other.z
  }                                                                                
                                                                                    
  pub fn add(&self, other: &Vector3) -> Vector3 {                              
    Vector3{                                                              
      x: self.x + other.x,
      y: self.y + other.y,
      z: self.z + other.z,
    }                                                                       
  }                                                                                
                                                                                    
  pub fn subtract(&self, other: &Vector3) -> Vector3 {                              
    Vector3{                                                              
      x: self.x - other.x,
      y: self.y - other.y,
      z: self.z - other.z,
    }                                                                       
  }         

  pub fn is_equal(&self, other: &Vector3, tol: f64) -> bool {
    (self.x - other.x).abs() < tol
    && (self.y - other.y).abs() < tol
    && (self.z - other.z).abs() < tol
  }

  pub fn cross(&self, other: &Vector3) -> Vector3{
    Vector3{                                                              
      x: self.y*other.z - self.z*other.y,
      y: self.z*other.x - self.x*other.z,
      z: self.x*other.y - self.y*other.x,
    }                                                                       
  }

  pub fn to_spherical_coordinates(&self) -> SphereVec3{

    let r = self.magnitude();
    let theta = (self.z/r).acos();
    let rho = r*f64::sin(theta);
    let mut phi = 0.0;
    
    if rho > 0.0{
      phi = (self.x/rho).acos();
     
      if self.y < 0.0 {
        phi = -phi;
      }
    }
    SphereVec3{
      r, theta, phi,
    }


  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
impl std::ops::Add<&Vector3> for &Vector3{
  type Output = Vector3;
  fn add(self, _rhs: &Vector3 ) -> Vector3 {
    self.add(_rhs)
  }
}

impl std::ops::Sub<&Vector3> for &Vector3{
  type Output = Vector3;
  fn sub(self, _rhs: &Vector3 ) -> Vector3 {
    self.subtract(_rhs)
  }
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[cfg(test)]
mod tests {
  use super::*;
  
  #[test]
  fn test_add(){
    let tol =1e-12; 
    let v = Vector3::from([1.0, 0.0, 1.0]);
    let u = Vector3::from([0.0, 1.0, 1.0]);
    let a = Vector3::from([1.0, 1.0, 2.0]);
    let b = &v + &u;
    assert!(a.is_equal(&b, tol));
  }

  #[test]
  fn test_sub(){
    let tol =1e-12; 
    let v = Vector3::from([1.0, 0.0, 1.0]);
    let u = Vector3::from([0.0, 1.0, 1.0]);
    let a = Vector3::from([1.0, -1.0, 0.0]);
    let b = &v - &u;
    assert!(a.is_equal(&b, tol));
  }
  
  #[test]
  fn test_magnitude(){
    let tol =1e-12; 
    let v = Vector3::from([1.0, 1.0, 1.0]);
    let a = v.magnitude();
    let b = f64::sqrt(3.0);

    assert!((a-b).abs() < tol);
  }

  #[test]
  fn test_normalize(){
    let tol =1e-12; 
    let v = Vector3::from([1.0, -10.0, 0.1]).normalize();
    let a = v.magnitude();
    let b = 1.0;
    
    assert!((a-b).abs() < tol);
  }
  
  #[test]
  fn test_scale(){
    let tol =1e-12; 
    let v = Vector3::from([1.0, -1.0, 2.0]);
    let u = v.scale(2.0);
    let a = Vector3::from([2.0, -2.0, 4.0]);
    assert!(a.is_equal(&u, tol));
  }

  
  #[test]
  fn test_dot(){
    let tol =1e-12; 
    let v = Vector3::from([1.0, -1.0, 1.0]);
    let u = Vector3::from([0.0, 1.0, 1.0]);
    let a = 0.0;
    let b = v.dot(&u);
    assert!((a-b).abs() < tol);
  }

  
  #[test]
  fn test_cross(){
    let tol =1e-12; 
    let x = Vector3::from([1.0, 0.0, 0.0]);
    let y = Vector3::from([0.0, 1.0, 0.0]);
    let z = Vector3::from([0.0, 0.0, 1.0]);

    assert!(z.is_equal(&x.cross(&y), tol));
    assert!(y.is_equal(&z.cross(&x), tol));
    assert!(x.is_equal(&y.cross(&z), tol));
  }

  #[test]
  fn test_to_spherical_coordinates(){
    let tol =1e-12; 
    let x = Vector3::from([1.0, 0.0, 0.0]);
    let y = Vector3::from([0.0, 1.0, 0.0]);
    let z = Vector3::from([0.0, 0.0, 1.0]);

    let rx = x.to_spherical_coordinates();
    let ry = y.to_spherical_coordinates();
    let rz = z.to_spherical_coordinates();

    assert!((rx.r-1.0).abs() < tol);
    assert!((rx.theta-PI/2.0).abs() < tol);
    assert!((rx.phi-0.0).abs() < tol);

    assert!((ry.r-1.0).abs() < tol);
    assert!((ry.theta-PI/2.0).abs() < tol);
    assert!((ry.phi-PI/2.0).abs() < tol);

    assert!((rz.r-1.0).abs() < tol);
    assert!((rz.theta-0.0).abs() < tol);
    assert!((rz.phi-0.0).abs() < tol);

  }


}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
