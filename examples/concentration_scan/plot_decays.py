import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

def main():

  dir_list = [
    "CluE-proton_frac_1",
    "CluE-proton_frac_0.75",
    "CluE-proton_frac_0.5",
    "CluE-proton_frac_0.25",
    "CluE-proton_frac_0",
  ]

  proton_fractions = [1.0, 0.75, 0.5, 0.25, 0.0]

  times = []
  signals = []
  for folder in dir_list:
    times.append( read_time_axis(f"{folder}/time_axis.csv") )
    signals.append( read_signal_file(f"{folder}/signal.csv") )
    

  plt.rc('font', size=12)
  plt.rc('axes', labelsize=12)

  fig, (ax) = plt.subplots(1,1,figsize=(8, 6) );
  cmap = matplotlib.colormaps['viridis']

  for ii, frac in enumerate(proton_fractions):
    ax.plot(times[ii]*1e6,np.real(signals[ii]),
      color = cmap(frac), linewidth=1);

  # Normalizer
  norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
  
  # creating ScalarMappable
  sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  sm.set_array([])

  plt.colorbar(sm, label=r"proton fraction (%)", orientation="vertical")  
  ax.set_xlabel(r"$2\tau$ (μs)");
  ax.set_ylabel("normalized echo amplitude");
  ax.set_ylim(0,1.0)
  ax.set_xscale("log")

  plt.savefig('CluE-sims.png')
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

  convert_to_python = {'i':'j', '\+-':'-'};
  v = (data[key].replace(convert_to_python,regex=True)).apply(
      lambda z: complex(z));
  return np.array(v)
#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------

if __name__ == "__main__":
  main();
