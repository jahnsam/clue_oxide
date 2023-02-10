use std::collections::HashSet;

pub fn are_vecs_equal<T: std::cmp::PartialEq>
(vec0: &[T], vec1: &[T])-> bool{
  if vec0.len() != vec1.len() {return false;}

  for ii in 0..vec0.len(){
    if vec0[ii] != vec1[ii]{ return false;}
  }
  true
}
//------------------------------------------------------------------------------
pub fn ceil(x: f64) -> f64{

  let mut a = x as i32;

  let err = (x-a as f64).abs();
  if err > 1e-12 && x >= 0.0{
    a += 1;
  }

  a as f64
}
//------------------------------------------------------------------------------

pub fn unique<T>(vec: Vec::<T>) -> Vec::<T>
where T: std::cmp::Eq  + std::hash::Hash
{
  vec.into_iter().collect::<HashSet<_>>()
      .into_iter()
      .collect()
}
//------------------------------------------------------------------------------


#[cfg(test)]
mod tests{
  use super::*;

  #[test]
  fn test_are_vecs_equal(){
    let a = vec![1,2,3];
    let b = vec![1,2,3];
    let c = vec![1,2,3,4];

    assert!(are_vecs_equal(&a,&b));
    assert!(!are_vecs_equal(&a,&c));
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_unique(){
    let a = vec![1,1,2,3,2];
    let mut a = unique(a);
    a.sort();
    assert_eq!(a,vec![1,2,3]);
  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_ceil(){
    assert_eq!(ceil(2.0),2.0);
    assert_eq!(ceil(2.1),3.0);
    assert_eq!(ceil(-2.1),-2.0);
    assert_eq!(ceil(-2.0),-2.0);
  }
}
