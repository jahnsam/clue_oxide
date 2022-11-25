use super::vector3::Vector3;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
pub struct Methyl{
  center: Vector3,
  normal: Vector3,
}

//------------------------------------------------------------------------------

impl Methyl{                                                                     
                                                                                 
pub fn from(r_carbon: Vector3,                                                   
    h0: Vector3,                                                                 
    h1: Vector3,                                                                 
    h2: Vector3,                                                                 
    ) -> Methyl{                                                                 
                                                                                 
  let center = (&( &h0 + &h1) + &h2).scale(1.0/3.0);                             
                                                                                 
  let x = (&h0 - &center).normalize();                                           
                                                                                 
  let y = (&h1 - &center).normalize();                                           
  let y = &y - &x.scale(x.dot(&y));                                              
                                                                                 
                                                                                 
  let mut normal = x.cross(&y);                                                  
                                                                                 
  let axis = &center -  &r_carbon;                                               
  if normal.dot(&axis) < 0.0 {                                                   
    normal = normal.scale(-1.0);                                                 
  }                                                                              
                                                                                 
  Methyl{ center, normal}                                                        
                                                                                 
}                                                                                
                                                                                 
pub fn center(&self) -> &Vector3 {                                                
  &self.center                                                           
}

pub fn normal(&self) -> &Vector3 {                                                
  &self.normal                                                      
}

}                                                                                
                                                                                    
                                                                                    
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


