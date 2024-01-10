import clue_oxide as clue
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

def main():
  options = clue.CluEOptions(
      config = {
        "input_structure_file": "\"MD_300K_25ns-30ns_center_0001.pdb\"",
        "radius": "25", # Å
        "detected_spin_position": "centroid_over_serials([28,29])",
        "detected_spin_g_matrix": "[2.0097, 2.0064, 2.0025]",
        "detected_spin_g_y": "diff(group(tempo_n) , group(tempo_o) )",
        "detected_spin_g_x": "diff(group(tempo_c1) , group(tempo_c19) )",
        "number_timepoints": "[101]",
        "time_increments": "[5e-8]", # s
        "cluster_method": "cce",
        "max_cluster_size": "3",
        "magnetic_field": "1.2", # T
        "pulse_sequence": "hahn",
        "neighbor_cutoff_3_spin_hahn_mod_depth": "3.23e-5",
        "neighbor_cutoff_3_spin_hahn_taylor_4": "1e17", # (rad/s)^4
      },
      group_list = [
        clue.Group("tempo_h", criteria = {
          "residues in": "[TEM]",
          "elements in": "[H]",
        }),
        clue.Group("tempo_n", criteria = {
          "residues in": "[TEM]",
          "elements in": "[N]",
        }),
        clue.Group("tempo_o", criteria = {
          "residues in": "[TEM]",
          "elements in": "[O]",
        }),
        clue.Group("tempo_c1", criteria = {
          "residues in": "[TEM]",
          "elements in": "[C]",
          "serials in": "[1]"
        }),
        clue.Group("tempo_c19", criteria = {
          "residues in": "[TEM]",
          "elements in": "[C]",
          "serials in": "[19]"
        }),
      ],
      spin_properties_list = [
        clue.SpinProperties("tempo_n", "14N", properties = {
          "hyperfine_coupling": "[20,20,100]*1e6", # Hz
          "hyperfine_x": "diff(bonded(tempo_c1),bonded(tempo_c19))",
          "hyperfine_y": "diff(particle,bonded(tempo_o))",
          "electric_quadrupole_coupling": "[-0.28, -1.47, 1.75]*1e6", # Hz
          "electric_quadrupole_x": "diff(bonded(tempo_c1),bonded(tempo_c19))",
          "electric_quadrupole_y": "diff(particle,bonded(tempo_o))",
        }),
      ])


  plt.rc('font', size=12)
  plt.rc('axes', labelsize=12)

  fig, (ax) = plt.subplots(1,1,figsize=(8, 6) );
#cmap = matplotlib.cm.get_cmap('viridis');
  cmap = matplotlib.colormaps['viridis']

  tunnel_splittings = np.array([0,20,40,80,100])*1e3;
  for nut in tunnel_splittings:
    options.set_spin_properties("tempo_h", "1H",
        key = "tunnel_splitting", value = f"{nut}")
    options.print()
    time , signal = options.run()
    ax.plot(np.array(time)*1e6,np.real(np.array(signal)),
      color = cmap(nut/max(tunnel_splittings)), linewidth=1);

  # Normalizer
  norm = matplotlib.colors.Normalize(vmin=0, vmax=max(tunnel_splittings)*1e-3)
  
  # creating ScalarMappable
  sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  sm.set_array([])

  plt.colorbar(sm, label=r"$\nu_\mathrm{t}$ (kHz)", orientation="vertical")  
  ax.set_xlabel(r"$2\tau$ (μs)");
  ax.set_ylabel("normalized echo amplitude");
  ax.axis([0,10,0,1.0] )

  plt.savefig('CluE-sims.png')
if __name__ == "__main__":
  main();
