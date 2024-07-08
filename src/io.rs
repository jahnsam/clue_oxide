use crate::CluEError;
use std::error::Error;
//------------------------------------------------------------------------------
/// This function writes signal data to a csv.
pub fn write_data<T>(data: &[Vec::<T>], filename: &str,
    headers: Vec::<String>) -> Result<(),CluEError> where T: std::fmt::Display
{
   if headers.len() != data.len(){
     return Err(CluEError::MissingHeader(filename.to_string()));
   }

   for datum in data.iter(){
     if datum.len() != data[0].len(){
       return Err(CluEError::AllVectorsNotSameLength(filename.to_string()));
     }
   }

    match write_vecs_to_csv(data,filename,headers){
      Ok(()) => Ok(()),
      Err(_) => Err(CluEError::CannotWriteFile(filename.to_string())),
    }
}
//------------------------------------------------------------------------------
// This function writes signal data to a csv.
// The function write_data calls this function, but has an extra error handling
// layer.
fn write_vecs_to_csv<T>(data: &[Vec::<T>],filename: &str,
    headers: Vec::<String>) 
  -> Result<(),Box<dyn Error>> where T: std::fmt::Display
{
  let n_data = data[0].len();

  let mut wtr = csv::Writer::from_path(filename)?;
  
  wtr.write_record(&headers)?;

  for ii in 0..n_data{
    
    let mut rec = Vec::<String>::with_capacity(n_data);
   
    for datum in data.iter(){
      rec.push(datum[ii].to_string());
    }

    wtr.write_record(&rec)?;
  }

  wtr.flush()?;
  
  Ok(())
}
//------------------------------------------------------------------------------
