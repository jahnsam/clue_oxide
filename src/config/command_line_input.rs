use crate::clue_errors::CluEError;
use substring::Substring;

/// The struct contains the processed command line input.
pub struct CommandLineInput{
  pub config_file: Option<String>,
  pub output_file: Option<String>,
  pub config_options: Vec::<String>,
  pub show_help: bool,
  pub show_license: bool,
  pub show_title: bool,
  pub show_version: bool,
  pub show_warrenty: bool,
  next_arg: NextArg
}

enum NextArg{
  InputConfig,
  Output,
  Unknown,
}

impl CommandLineInput{

  pub fn new(args: Vec::<String>) -> Result<Self,CluEError> {
    let mut cli = CommandLineInput{
      config_file: None,
      output_file: None,
      config_options: Vec::<String>::new(),
      show_help: false,
      show_license: false,
      show_title: true,
      show_version: false,
      show_warrenty: false,
      next_arg: NextArg::InputConfig,
    };

    let mut it = args.iter();
    it.next();
    for arg in it{

      if (*arg).substring(0,2)== "--"{
        cli.parse_long_option(arg)?;
        continue;
      }
      
      if(*arg).substring(0,1)== "-"{
        cli.parse_short_options((*arg).substring(1,(*arg).len()))?;
        continue;
      }

      cli.parse_next_argument(arg)?;
    }
    Ok(cli)
  }
  //----------------------------------------------------------------------------
  fn parse_next_argument(&mut self, arg: &str) -> Result<(),CluEError>{
    match self.next_arg{
      NextArg::InputConfig => self.config_file = Some(arg.to_string()),
      NextArg::Output => self.output_file = Some(arg.to_string()),
      NextArg::Unknown =>
        return Err(CluEError::UnrecognizedOption(arg.to_string() )),
    }

    if self.config_file.is_none() {
      self.next_arg = NextArg::InputConfig;
    }else if self.output_file.is_none(){
      self.next_arg = NextArg::Output;
    }else{
      self.next_arg = NextArg::Unknown;
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  fn parse_long_option(&mut self, option: &str) -> Result<(),CluEError>{
  
    match option {
      "--help" => self.show_help = true,
      "--license" => self.show_license = true,
      "--hide-title" => self.show_title = false,
      "--output" => self.next_arg = NextArg::Output,
      "--version" => self.show_version = true,
      "--warrenty" => self.show_warrenty = true,
      _ => return Err(CluEError::UnrecognizedOption(option.to_string() )),  
    }
    Ok(())
  }
  //----------------------------------------------------------------------------
  fn parse_short_options(&mut self, options: &str) ->  Result<(),CluEError>{
    
    let n_opts = options.len();

    for ii in 0..n_opts{
      let opt = options.substring(ii,ii+1);

      self.parse_short_option(opt)?;
    }
    Ok(())
  
  }
  //----------------------------------------------------------------------------
  fn parse_short_option(&mut self, option: &str) ->  Result<(),CluEError>{

    match option{
      "h" => self.show_help = true,
      "l" => self.show_license = true,
      "H" => self.show_title = false,
      "o" => self.next_arg = NextArg::Output,
      "V" => self.show_version = true,
      "W" => self.show_warrenty = true,
      _ => return Err(CluEError::UnrecognizedOption(option.to_string() )),  
    }
    Ok(())
  }  
  //----------------------------------------------------------------------------
}
