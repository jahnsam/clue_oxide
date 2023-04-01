from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

color_palette = [
  [222/255,66/255,91/255,1],
  [5/255,135/255,0, 1],
  [84/255,93/255,187/255,1],
  [255/255,166/255,0/255,1],
  [221/255,119/255,136/255,1],
  [102/255,153/255,204/255,1],
  [144/255,238/255,144/255,1],
  [253/255,230/255,47/255,1],
]
#-------------------------------------------------------------------------------
def main():
  if len(sys.argv) == 1:
    print("Please supply a csv file:")
    print("  python3 signal.csv")
    return

  signal = read_signal_file(sys.argv[1])
  time = np.linspace(0,1,len(signal))

  fig, (ax) = plt.subplots(1,1,figsize=(8, 6) );
  ax.plot(time,np.real(signal),
      linewidth=1,color=[5/255,135/255,0]);

  plt.show()

#-------------------------------------------------------------------------------
def read_signal_file(csv_file):
  data = pd.read_csv(csv_file);
  #t = data['time'];
  convert_to_python = {'i':'j', '\+-':'-'};
  v = (data['signal'].replace(convert_to_python,regex=True)).apply(
      lambda z: complex(z));
  return np.array(v)
#-------------------------------------------------------------------------------


#===============================================================================
if __name__ == "__main__":
  main()
#===============================================================================
