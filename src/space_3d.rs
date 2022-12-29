
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug, Clone)]
pub struct SymmetricTensor3D{
  elements: [f64; 6],
}

impl SymmetricTensor3D{
  pub fn zeros() -> Self{
    SymmetricTensor3D{elements: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}
  }

  pub fn xx(&self) -> f64 { self.elements[0]}
  pub fn xy(&self) -> f64 { self.elements[1]}
  pub fn xz(&self) -> f64 { self.elements[2]}
  pub fn yx(&self) -> f64 { self.xy() }
  pub fn yy(&self) -> f64 { self.elements[3]}
  pub fn yz(&self) -> f64 { self.elements[4]}
  pub fn zx(&self) -> f64 { self.xz()}
  pub fn zy(&self) -> f64 { self.yz()}
  pub fn zz(&self) -> f64 { self.elements[5]}

  pub fn set_xx(&mut self, value: f64) { self.elements[0] = value;}
  pub fn set_xy(&mut self, value: f64) { self.elements[1] = value;}
  pub fn set_xz(&mut self, value: f64) { self.elements[2] = value;}
  pub fn set_yx(&mut self, value: f64) { self.elements[1] = value;}
  pub fn set_yy(&mut self, value: f64) { self.elements[3] = value;}
  pub fn set_yz(&mut self, value: f64) { self.elements[4] = value;}
  pub fn set_zx(&mut self, value: f64) { self.elements[2] = value;}
  pub fn set_zy(&mut self, value: f64) { self.elements[4] = value;}
  pub fn set_zz(&mut self, value: f64) { self.elements[5] = value;}

}

impl ToString for SymmetricTensor3D{
  fn to_string(&self) -> String{
    format!("[{}, {}, {}, {}, {}, {}]", self.elements[0], self.elements[1],
        self.elements[2], self.elements[3],self.elements[4], self.elements[5])
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#[derive(Debug, Clone)]
pub struct Vector3D{
  elements: [f64; 3],
}

impl Vector3D{
  pub fn zeros() -> Self{
    Vector3D{elements: [0.0, 0.0, 0.0]}
  }

  pub fn x(&self) -> f64 { self.elements[0]}
  pub fn y(&self) -> f64 { self.elements[1]}
  pub fn z(&self) -> f64 { self.elements[2]}

  pub fn set_x(&mut self, value: f64) { self.elements[0] = value;}
  pub fn set_y(&mut self, value: f64) { self.elements[0] = value;}
  pub fn set_z(&mut self, value: f64) { self.elements[0] = value;}
}

impl ToString for Vector3D{
  fn to_string(&self) -> String{
    format!("[{}, {}, {}]", self.elements[0], self.elements[1],
        self.elements[2])
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
