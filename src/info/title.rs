use crate::info::version;

pub fn print_title() {  
println!("\
##############################################################################\n\
#                                                                            #\n\
#                      ______    _             ________                      #\n\
#                     /  ____|  | |           | _______|                     #\n\
#                    |  |       | |   _   _   | |_____                       #\n\
#                    |  |       | |  | | | |  |  _____|                      #\n\
#                    |  |____   | |  | |_| |  | |______                      #\n\
#                     \\______|  |_|  \\_____/  |________|                     #\n\
#                                                                            #\n\
#                 ___               ______               ___                 #\n\
#                    \\             /      \\             /                    #\n\
#                     \\           /        \\           /                     #\n\
#                      \\    O    /          \\    O    /                      #\n\
#                       \\       /            \\       /                       #\n\
#                        \\     /              \\     /                        #\n\
#                ()       \\   /       ()       \\   /       ()                #\n\
#                          \\ /                  \\ /                          #\n\
#                           V                    V                           #\n\
#                                                                            #\n\
##############################################################################\n");
  println!("CLuE Copyright (C) 2022 Samuel M. Jahn");
  version::print_version();
  println!("\nThis program comes with ABSOLUTELY NO WARRANTY; for details run \
\"clue --warrenty\".");
  println!("This is free software, and you are welcome to redistribute it");
  println!("under certain conditions; run \"clue --license\" for details.\n");
  println!("For usage information; run \"clue --help\" for more details.\n");


}
