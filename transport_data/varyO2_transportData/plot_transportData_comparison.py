#Python script to visualize the various transport data(mobility, ionization coeff., attachment coeff., etc when varying the water vapour percentage)
# for two different transport data files
import matplotlib.pyplot as plt
import numpy as np
import re


def get_transportdata(fname, dataname):
    with open(fname, "r") as f:
        fr = f.read()
    fr = re.sub(r'\r\n', "\n", fr)
    lines = fr.splitlines()
    ix = lines.index(r"{}".format(dataname))+3 
    i1 = ix
    while lines[i1].strip() != "":
        i1 += 1
    tbl_string = ";".join(lines[ix:i1-1])
    return np.matrix(tbl_string)



#Computing the gas density to unscale the coeffcients
p = 1e5 #Pascal
T = 300 #Kelvin
kb = 1.38e-23 #m^2 kg s^-2 K^-1
N = p/(kb*T)
old = "air_chemistry_v2_old.txt"
new = "air_chemistry_v2_new.txt"


coeff = ["Mobility *N (1/m/V/s)", "Diffusion coefficient *N (1/m/s)", "Townsend ioniz. coef. alpha/N (m2)", "Townsend attach. coef. eta/N (m2)"]
coeff_label = ["Mobility (m2/V/s)", "Diffusion coefficient (m2/s)", "Townsend ioniz. coef. alpha (1/m)", "Townsend attach. coef. eta (1/m)"]
scale = [N**-1, N**-1, N, N, 1.0]

for i1, c in enumerate(coeff):
    plt.figure(i1)
    of = get_transportdata(old, c)
    nf = get_transportdata(new, c)
    plt.plot(of[:,0], of[:,1]*scale[i1], label="old "+coeff_label[i1])
    plt.plot(nf[:,0], nf[:,1]*scale[i1], label="new "+coeff_label[i1])
    plt.title(coeff_label[i1])
    plt.legend()
    plt.xlabel("E/N (Td)")
    #plt.xscale("log")

plt.show()


