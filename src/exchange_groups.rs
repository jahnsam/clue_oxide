use super::vector3::Vector3;

pub trait Translate{ fn translate(&mut self, r: &Vector3); }

pub trait GetIndices{fn indices(&self) -> Vec::<usize>; }

pub trait GetExchangeCoupling{fn exchange_coupling(&self) -> f64; }


#[derive(Debug, Clone)]
pub enum ExchangeGroup{
  Methyl(C3Rotor),
  PrimaryAmonium(C3Rotor),
}

impl GetIndices for ExchangeGroup{
  fn indices(&self) -> Vec::<usize>{
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.indices(),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.indices(),
    }
  }
}

impl Translate for ExchangeGroup{
  fn translate(&mut self, r: &Vector3){
    match self{
      ExchangeGroup::Methyl(rotor) => rotor.translate(r),
      ExchangeGroup::PrimaryAmonium(rotor) => rotor.translate(r),
    }
  }
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#[derive(Debug, Clone)]
pub struct C3Rotor{
  center: Vector3,
  normal: Vector3,
  indices: [usize;3],
}

//------------------------------------------------------------------------------

impl C3Rotor{                                                                     
                                                                                 
pub fn from(r_carbon: Vector3,                                                   
    h0: Vector3,                                                                 
    h1: Vector3,                                                                 
    h2: Vector3,
    indices: [usize;3],    
    ) -> Self{                                                                 
                                                                                 
  let center = (&( &h0 + &h1) + &h2).scale(1.0/3.0);                             
                                                                                 
  let x = (&h0 - &center).normalize();                                           
                                                                                 
  let y = (&h1 - &center).normalize();                                           
  let y = &y - &x.scale(x.dot(&y));                                              
                                                                                 
                                                                                 
  let mut normal = x.cross(&y);                                                  
                                                                                 
  let axis = &center -  &r_carbon;                                               
  if normal.dot(&axis) < 0.0 {                                                   
    normal = normal.scale(-1.0);                                                 
  }                                                                              
                                                                                 
  C3Rotor{ center, normal, indices}                                                        
                                                                                 
}                                                                                
                                                                                 
pub fn center(&self) -> &Vector3 {                                                
  &self.center                                                           
}

pub fn normal(&self) -> &Vector3 {                                                
  &self.normal                                                      
}
}

impl GetIndices for C3Rotor{
  fn indices(&self) -> Vec::<usize>{
    
    let mut out = Vec::<usize>::with_capacity(3);

    for h in self.indices{
      out.push(h);
    }
    out
  }
}

impl Translate for C3Rotor{
  fn translate(&mut self, r: &Vector3){
   self.center = &self.center + r;
  }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

