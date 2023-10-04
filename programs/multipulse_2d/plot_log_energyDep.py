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

pref = "./multiPulse/varyTinter/bench_27kV_"
suf = "mus_log.txt"
log_files = [pref+x+suf for x in  ["10", "25", "50", "200", "500"]]
#log_files.insert(3, "./multiPulse/bench_27kV_eta0p3_log.txt")
logs = [pd.read_csv(f, delim_whitespace=True) for f in log_files]
log_labels = ["100kHz", "40kHz", "20kHz", "5kHz", "2kHz"]
#log_pulses = [20, 20, 20, 13, 7, 3]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
#fig.suptitle('\n'.join(numbered_files))

print(logs[0].columns)
for i, log in enumerate(logs):
    for y in args.y:
        # a = args.yscale*np.diff(log[y])
        # print(a[a > 1.0])
        #print(log_files[i])
        #a = args.yscale*np.max(log[y])/ log_pulses[i]
        #print(args.yscale*np.max(log[y])/ log_pulses[i])
        axes.plot(args.xscale*log[args.x], args.yscale*log[y], label=log_labels[i])

plt.xlabel(args.x+" (milli sec.)")
plt.ylabel("Maximum temperature (K)")
plt.legend()
plt.show()
