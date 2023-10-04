#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("-x", type=str, default='time', help="Name of x variable")
p.add_argument("-y", type=str, nargs='+', default='[max(E)]',
               help="Name of y variables")
p.add_argument("-yscale", type=float, default=1.0, help="value to multiply the y-variable")
p.add_argument("-xscale" ,type=float, default=1.0, help="value to multiply the x-variable")
args = p.parse_args()

pref = "./multiPulse/bench_"
suf = "kV_log.txt"
log_files = [pref+x+suf for x in  ["22", "24", "27", "32"]]
log_files.append("./multiPulse/bench_30kV_20_pulses_log.txt")
logs = [pd.read_csv(f, delim_whitespace=True) for f in log_files]
log_labels = ["22kV", "24kV", "27kV", "32kV", "30kV"]
log_pulses = [13, 13, 13, 9, 12]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
#fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    for y in args.y:
        # a = args.yscale*np.diff(log[y])
        # print(a[a > 1.0])
        a = args.yscale*np.max(log[y])/ log_pulses[i]
        print(args.yscale*np.max(log[y])/ log_pulses[i])
        axes.plot(args.xscale*log[args.x], args.yscale*log[y], label=log_labels[i])

plt.xlabel(args.x+" (milli sec.)")
plt.ylabel("Total energy deposited ($\mu$ J)")
plt.legend()
plt.show()
