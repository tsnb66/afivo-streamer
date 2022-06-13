#Script to generate a set of scaled input values for simulations in afivo-streamer. In this particular script, we input an array of pressures, and the output is a set of scaled variables that need to be updated in the .cfg file
import numpy as np



#Reference simulation parameters
V_ref = 20e3 #Voltage in V
P_ref = 80e-3
d_ref = 160e-3 #Domain gap length in m
BG_ref = 1e11 #Background density in m^-3
ref_electrode_dx_ref = 2e-5
rod_radius_ref = 4e-4
tip_radius_ref = 5e-5
refine_maxdx_ref = 4e-3
refine_mindx_ref = 1e-7
#Input pressure array (in bar)
P_out = [133e-3, 233e-3, 400e-3]



#Scaled domain gaps, so that E/N is constant, and we assume V is constant

d_out = [(d_ref*P_ref)/x for x in P_out]


#Scaling the background density (ONLY FOR 2D CASES)

BG_out = [BG_ref*(x/P_ref)**2 for x in P_out]


#Scaling the electrode and the dx values around the electrode

rod_rad_out = [rod_radius_ref*P_ref/x for x in P_out]
tip_rad_out = [tip_radius_ref*P_ref/x for x in P_out]
ref_el_dx_out = [ref_electrode_dx_ref*P_ref/x for x in P_out]

#Rescaling the dx control values
refine_maxdx_out = [refine_maxdx_ref*P_ref/x for x in P_out]
refine_mindx_out = [refine_mindx_ref*P_ref/x for x in P_out]


# Printing the values
print("P out: ", P_out)
print("D: ", d_out)
print("BG: ", BG_out)
print("Rod Radius: ", rod_rad_out)
print("Tip Radius: ", tip_rad_out)
print("Electrode refinement dx: ", ref_el_dx_out)
print(" Max dx: ", refine_maxdx_out)
print(" Min dx: ", refine_mindx_out)


