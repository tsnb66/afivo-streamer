#!/usr/bin/env python3

# Author: Jannis Teunissen

import argparse
import re
import numpy as np
import sys


def getArgs():
    "Set the arguments, read them in, and return them"
    parser = argparse.ArgumentParser(description="Converter for Bolsig+ data",
                formatter_class= argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", type=str, help="input file")
    parser.add_argument("outfile", type=str, help="output file")
    parser.add_argument("-transport", type=str, nargs='+',
                        default=["Mean energy (eV)",
                                 "Mobility *N (1/m/V/s)",
                                 "Diffusion coefficient *N (1/m/s)",
                                 "Townsend ioniz. coef. alpha/N (m2)",
                                 "Elastic power loss /N (eV m3/s)",
                                 "Inelastic power loss /N (eV m3/s)",
                                 "Townsend attach. coef. eta/N (m2)"],
                        help="Which transport coefficients to include")
    #parser.add_argument("-transport", type=str, nargs='+',
    #                    default=["Mean energy (eV)",
    #                             "Mobility *N (1/m/V/s)",
    #                             "Diffusion coefficient *N (1/m/s)",
    #                             "Townsend ioniz. coef. alpha/N (m2)",
    #                             "Townsend attach. coef. eta/N (m2)"],
    #                    help="Which transport coefficients to include")
    parser.add_argument("-colltypes", type=str, nargs='+',
                        default=["Ionization", "Attachment"],
                        help="Which collision types to include")
    parser.add_argument("-colls", type=str, nargs='+', default=[],
                        help="Other collisions to include (e.g. C2 C5)")
    return parser.parse_args()


def convert():
    cfg = getArgs()

    with open(cfg.infile, 'r') as f:
        fr = f.read()
    fr = re.sub(r'\r\n', '\n', fr)  # Fix newlines
    lines = fr.splitlines()

    transport_labels = ['R#', 'E/N']
    ix = lines.index(r' Transport coefficients') + 1
    while lines[ix].startswith(r'A'):
        tmp = lines[ix].split()
        transport_labels.append(' '.join(tmp[1:]))
        ix += 1

    # Skip one line with column labels
    i0 = ix + 1

    # Find end of table
    i1 = i0
    while lines[i1].strip() != "":
        i1 += 1

    tbl_string = ";".join(lines[i0:i1])
    transport_data = np.matrix(tbl_string)

    rate_labels = ['R#', 'E/N', 'eV']
    rate_desc = [None, None, None]
    ix = lines.index("Rate coefficients (m3/s)") + 1
    while lines[ix].startswith(r' C'):
        tmp = lines[ix].split()
        rate_labels.append(tmp[0])
        rate_desc.append(' '.join(tmp))
        if tmp[2] in cfg.colltypes:
            cfg.colls.append(tmp[0])
        ix += 1

    # Skip one line with column labels
    i0 = ix + 1

    # Find end of table
    i1 = i0
    while lines[i1].strip() != "":
        i1 += 1

    tbl_string = ";".join(lines[i0:i1])
    rate_data = np.matrix(tbl_string)
    #--------------------------Energy loss coeffs-----------------
#    energyLoss_labels = ['R#', 'E/N', 'eV']
#    energyLoss_desc = [None, None, None]
#    ix = lines.index("Energy loss coefficients (eV m3/s))") + 1
#    while lines[ix].startswith(r' C'):
#        tmp = lines[ix].split()
#        energyLoss_labels.append(tmp[0])
#        energyLoss_desc.append(' '.join(tmp))
#        if tmp[2] in cfg.colltypes:
#            cfg.colls.append(tmp[0])
#        ix += 1

#    # Skip one line with column labels
#    i0 = ix + 1

#    # Find end of table
#    i1 = i0
#    while lines[i1].strip() != "":
#        i1 += 1

#    tbl_string = ";".join(lines[i0:i1])
#    energyLoss_data = np.matrix(tbl_string)
    #-----------------------------------

    # Now write output
    with open(cfg.outfile, 'wb') as f:
        ix = transport_labels.index("E/N")
        for name in cfg.transport:
            if name in transport_labels:
                iy = transport_labels.index(name)
                write_entry(name, transport_data[:, ix],
                            transport_data[:, iy], f)
            else:
                print("Warning: {} not found".format(name))
                write_entry(name, transport_data[:, ix],
                            np.zeros(transport_data.shape[0]), f,
                            "Warning: no data was found")

        ix = rate_labels.index("E/N")
        for name in cfg.colls:
            if name in rate_labels:
                iy = rate_labels.index(name)
                write_entry(rate_desc[iy], rate_data[:, ix],
                            rate_data[:, iy], f)
            else:
                print("Warning: {} not found".format(name))

#        ix = energyLoss_labels.index("E/N")
#        for name in cfg.colls:
#            if name in energyLoss_labels:
#                iy = energyLoss_labels.index(name)
#                write_entry(energyLoss_desc[iy], energyLoss_data[:, ix],
#                            energyLoss_data[:, iy], f)
#            else:
#                print("Warning: {} not found".format(name))                



def write_entry(entry_name, x_data, y_data, f, comment=None):
    out_data = np.column_stack([x_data, y_data])
    
    print("Writing: ", entry_name)
    hdr = '\n\n' + entry_name + '\n' + 'COMMENT: generated by bolsig_convert.py'
    
    if comment:
    	hdr += "\nCOMMENT: " + comment
    	
    hdr += '\n-----------------------\n'
    
    
    ftr = '-----------------------\n'
    f.write(hdr.encode('ascii'))
    np.savetxt(f, out_data)
    f.write(ftr.encode('ascii'))


if __name__ == '__main__':
    convert()
