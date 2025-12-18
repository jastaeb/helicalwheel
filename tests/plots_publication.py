
#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path  
import os




#%%
import helicalwheel as hw

# %% Inputting the sequences  
aap2_seq_wt = "LVWELIRWLQAVAHQWQTIT"
aap2_seq_pm = "AAAEAARAAQAAAHQAQTAT" 


# Check that sequences are of the same length
if len(aap2_seq_wt) != len(aap2_seq_pm):
    raise ValueError("Wild-type and mutant sequences must be of the same length.")
else:
    mut_res = [i != y for i, y in zip(aap2_seq_wt, aap2_seq_pm)]

# %% Creating dataframes with helix coordinates
df_wt=hw.helix_coordinates(aap2_seq_wt)
df_wt_rotated=hw.helix_coordinates(aap2_seq_wt, rotate_deg=180)
df_pm=hw.helix_coordinates(aap2_seq_pm)
df_pm_rotated=hw.helix_coordinates(aap2_seq_pm, rotate_deg=180)

# Adding mutation information to the dataframes
df_wt["mutated"]=mut_res
df_wt_rotated["mutated"]=mut_res
df_pm["mutated"]=mut_res
df_pm_rotated["mutated"]=mut_res

#%% Plotting the helical wheels
first_pos = 21 # position in total protein sequence of first helical residue 

saving_path = Path("/Users/jann/Documents/GitHub/helicalwheel/plots")
saving_path.mkdir(parents=True, exist_ok=True)

ft1=0
lt1=10
ft2=10
lt2=len(aap2_seq_wt)


# 1. Wild-type AAP2
# Entire helix
wt_plot_1 = hw.helixplot(df_wt, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_wt.svg"), format="svg", bbox_inches='tight', pad_inches=0)


# Two individual plots figure assembly
wt_plot_1 = hw.helixplot(df_wt.iloc[ft1:lt1], legend=False, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_wt_{ft1+first_pos}_{lt1+first_pos}.svg"), format="svg", bbox_inches='tight', pad_inches=0)

wt_plot_1_rotated = hw.helixplot(df_wt_rotated.iloc[ft1:lt1], legend=False, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_wt_{ft1+first_pos}_{lt1+first_pos}_rotated.svg"), format="svg", bbox_inches='tight', pad_inches=0)

wt_plot_2 = hw.helixplot(df_wt.iloc[ft2:], legend=False, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_wt_{ft2+first_pos}_{lt2+first_pos}.svg"), format="svg", bbox_inches='tight', pad_inches=0)

wt_plot_2_rotated = hw.helixplot(df_wt_rotated.iloc[ft2:], legend=False, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_wt_{ft2+first_pos}_{lt2+first_pos}_rotated.svg"), format="svg", bbox_inches='tight', pad_inches=0)


# 2. Mutant AAP2
# Entire helix
pm_plot_1 = hw.helixplot(df_pm, first_pos=first_pos, makebold=True)
# Two individual plots figure assembly
pm_plot_1 = hw.helixplot(df_pm.iloc[0:10], legend=False, makebold=True, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_pm_{ft1+first_pos}_{lt1+first_pos}.svg"), format="svg", bbox_inches='tight', pad_inches=0)

pm_plot_1_rotated = hw.helixplot(df_pm_rotated.iloc[0:10], legend=False, makebold=True, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_pm_{ft1+first_pos}_{lt1+first_pos}_rotated.svg"), format="svg", bbox_inches='tight', pad_inches=0)

pm_plot_2 = hw.helixplot(df_pm.iloc[10:], legend=False, makebold=True, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_pm_{ft2+first_pos}_{lt2+first_pos}.svg"), format="svg", bbox_inches='tight', pad_inches=0)

pm_plot_2_rotated = hw.helixplot(df_pm_rotated.iloc[10:], legend=False, makebold=True, first_pos=first_pos)
plt.savefig(os.path.join(saving_path, f"helical_wheel_pm_{ft2+first_pos}_{lt2+first_pos}_rotated.svg"), format="svg", bbox_inches='tight', pad_inches=0)

# %%
