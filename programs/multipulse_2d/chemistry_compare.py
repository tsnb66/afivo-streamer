#!/usr/bin/env python3

# Authors: Hemaditya Malla, Jannis Teunissen

import numpy as np
import matplotlib.pyplot as plt

def plot_soistuff(fname, time_interval, soi, symb):
    # Assume the other files are in the same folder
    base_name = fname.replace('_rates.txt', '')

    with open(base_name + '_species.txt', 'r') as f:
        species_list = [x.strip() for x in f.readlines() if x.strip()]
    with open(base_name + '_reactions.txt', 'r') as f:
        reactions_list = [x.strip() for x in f.readlines() if x.strip()]

    stoich_matrix = np.loadtxt(base_name + '_stoich_matrix.txt')
    n_species, n_reactions = stoich_matrix.shape

    tmp = np.loadtxt(base_name + '_rates.txt')
    time = tmp[:, 0]
    rates = tmp[:, 1:]

    # Only consider chemistry within the given interval
    t1_idx = np.where(time >= time_interval[0])[0][0]
    # +1 for array slicing, but do not go beyond array size
    t2_idx = np.where(time <= time_interval[1])[0][-1] + 1
    t2_idx = min(t2_idx, len(time))

    time = time[t1_idx:t2_idx]
    time = (time-time_interval[0])/(time_interval[1]-time_interval[0])
    rates = rates[t1_idx:t2_idx]

    # Subtract initial state from rates
    rates = rates - rates[0]
    # Visualize the source and sink reactions for a given specie
    sidx = species_list.index(soi)

    srce_idx = np.where(stoich_matrix[sidx, :] > 0)[0]
    sink_idx = np.where(stoich_matrix[sidx, :] < 0)[0]
    titles = ['Source', 'Sink']

    for i, (ix, text) in enumerate(zip([srce_idx, sink_idx], titles)):
        amount = stoich_matrix[sidx, ix] * rates[:, ix]
        frac = amount[-1]/amount[-1].sum()

        threshold = 0.01
        for j, idx in enumerate(ix):
            if frac[j] > threshold:
                ax[i].plot(time, amount[:, j], label=reactions_list[idx] +
                           f' ({100*frac[j]:.2f}%)', marker=symb)

        ax[i].set_title(text + ' reactions')
        ax[i].set_xlabel('Time (s)')
        ax[i].set_ylabel('Production (#)')
        #ax[i].legend(bbox_to_anchor=(1.1, 0.5))
        ax[i].legend()

    gross_prod = np.dot(rates[:, srce_idx], stoich_matrix[sidx, srce_idx])
    net_prod = np.dot(rates, stoich_matrix[sidx])
    ax[2].plot(time, gross_prod, label='gross production', marker=symb)
    ax[2].plot(time, net_prod, label='net production', marker=symb)
    ax[2].set_xlabel('Time (s)')
    ax[2].set_ylabel('Production (#)')
    ax[2].legend()
    fig.suptitle(f'{len(srce_idx)+len(sink_idx)} of {n_reactions}'
                 f' influence {soi}')
    return True




fig, ax = plt.subplots(3, figsize=(8, 12), sharex=True)
#p1 = plot_soistuff("newBGChem/20O2/20O2_0p5mus_rates.txt", [620e-9,640e-9], "e", "o")
#tinterval = [81e-9, 1083e-9]
#p1 = plot_soistuff("varyRiseFall_133mbar/20O2_0p5mus_1nsrf_rates.txt", [582e-9, 584e-9], "e", ".")
#p1 = plot_soistuff("varyRiseFall_133mbar/20O2_0p5mus_10nsrf_rates.txt", [600e-9, 610e-9], "e", "^")
#tinterval = [81e-9, 1083e-9]
#p1 = plot_soistuff("varyRiseFall_133mbar/20O2_1mus_1nsrf_rates.txt", tinterval, "e", "^")
tinterval = [40e-9, 55e-9]
#p1 = plot_soistuff("1bar_1ns_newerChem/20O2/20O2_500ns_rates.txt", tinterval, "e", "o")
#p1 = plot_soistuff("output/20O2_500ns_retest_rates.txt", tinterval, "e", "^")
p1 = plot_soistuff("1bar_1ns_newerChem/20O2/20O2_25ns_rates.txt", tinterval, "e", "o")
p1 = plot_soistuff("output/20O2_25ns_retest_master_rates.txt", tinterval, "e", "^")
#p1 = plot_soistuff("output/20O2_500ns_retestOldChem_rates.txt", tinterval, "e", "*")
#p1 = plot_soistuff("newerBGChem/20O2_1mus_1nsrf_rates.txt", [81e-9, 1083e-9], "e", "o")
#p1 = plot_soistuff("newerBGChem/20O2_1mus_10nsrf_rates.txt", [90e-9, 1110e-9], "e", "^")
#p1 = plot_soistuff("newerBGChem/20O2_1mus_20nsrf_rates.txt", [100e-9, 1140e-9], "e", "*")
plt.tight_layout()
plt.show()
