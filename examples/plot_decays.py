import matplotlib
matplotlib.use('TkAgg')
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
  if len(sys.argv) < 3:
    print("Please supply a csv file:")
    print("  python3 time_axis.csv signal.csv")
    return

  time = read_time_axis(sys.argv[1])
  signals = []

  for ii in range(2,len(sys.argv)):
    signal = read_signal_file(sys.argv[ii])
    signals.append(signal)

  fig, (ax) = plt.subplots(1,1,figsize=(8, 6) );
  
  for ii,signal in enumerate(signals):
    color = color_palette[ii % 8]
    ax.plot(time*1e6,np.real(signal),
        linewidth=1,color=color);

  ax.set_xlabel("time (μs)")
  ax.set_ylabel("Re(signal)")
  plt.show()

#-------------------------------------------------------------------------------
def read_time_axis(csv_file):
  data = pd.read_csv(csv_file);
  t = data['time_axis'];
  return np.array(t);
#-------------------------------------------------------------------------------
def read_signal_file(csv_file):
  data = pd.read_csv(csv_file);
  order = 1
  key = "signal"
  while True:
    new_key = f"signal_{order}"
    if new_key in data:
      key = new_key
      order += 1
    else:
      break;

  #t = data['time'];
  convert_to_python = {'i':'j', '\+-':'-'};
  v = (data[key].replace(convert_to_python,regex=True)).apply(
      lambda z: complex(z));
  return np.array(v)
#-------------------------------------------------------------------------------


#===============================================================================
if __name__ == "__main__":
  main()
#===============================================================================
