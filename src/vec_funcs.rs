use std::collections::HashSet;


pub fn unique<T>(vec: Vec::<T>) -> Vec::<T>
where T: std::cmp::Eq  + std::hash::Hash
{
  vec.into_iter().collect::<HashSet<_>>()
      .into_iter()
      .collect()
}
