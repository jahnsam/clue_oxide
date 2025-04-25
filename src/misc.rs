
//------------------------------------------------------------------------------
pub fn eq_variant<T>(a: &T, b: &T) -> bool
{
  std::mem::discriminant(a) == std::mem::discriminant(b)
}
//------------------------------------------------------------------------------
pub fn are_all_same_type<T>(array: &Vec::<T>) 
  -> bool
{
  if array.len() <= 1 { return true; }

  for value in array.iter(){
    if !eq_variant(value,&array[0]){
      return false;
    } 
  }
  true
}
//------------------------------------------------------------------------------

