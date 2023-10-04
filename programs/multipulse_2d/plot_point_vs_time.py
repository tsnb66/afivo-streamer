
#!/usr/bin/env python3

import numpy as np
import argparse
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("pointVtime", type=str, nargs='+', help="Input point_vs_time file(s)")
#p.add_argument("-x", type=str, default='time', help="Name of x variable")
#p.add_argument("-y", type=str, default='time', default='[max(E)]',
#               help="Name of y variables")
p.add_argument("-yscale", type=float, default=1.0, help="value to multiply the y-variable")
p.add_argument("-xscale" ,type=float, default=1.0, help="value to multiply the x-variable")
args = p.parse_args()

logs = [np.loadtxt(f) for f in args.pointVtime]
numbered_files = [f'{i}: {f}' for i, f in enumerate(args.pointVtime)]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    axes.plot(args.xscale*log[:,0], args.yscale*log[:,1], label=f"{i}")

plt.xlabel("time")
plt.legend()
plt.show()
