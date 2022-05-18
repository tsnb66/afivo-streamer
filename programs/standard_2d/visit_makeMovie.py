#!/usr/bin/env python3

# This script saves a png for a given variable
#  for a given time iterations

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys
import numpy as np

def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to plot a surface plot in png for a given variable''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog='visit -nowin -cli -s visit_savefigure.py')
    pr.add_argument('database', type=str,
                    help='Database name (e.g. sim.silo or "sim_*.silo"')
    pr.add_argument('variable', type=str,
                    help='Name of variable')
    pr.add_argument('-scale', type=str,
                    help='The type of scaling to use (Linear, Log, Skewed)')
    pr.add_argument('-pltLims', nargs="+", type=float, 
    default=[float("NaN"), float("NaN")], 
                    help='Lower and max. limits of the scale. Uses the max, min value of the corresponding plot if not specified, wont save the figure if you have a non-positive val for Log scaling')
    pr.add_argument('-i0', type=int,
                    help='Start index in database (0 to N-1 for N files)')
    pr.add_argument('-i1', type=int,
                    help='Stop index in database (0 to N-1 for N files)')
    pr.add_argument('-outputDir', type=str, default='.',
                    help='Base path name for output directory')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    try: 
        import visit as v
    except ImportError:
        pr.print_help()
        raise

    args = pr.parse_args()
    

    

    # Open either a database or a single file
    use_database = ('*' in args.database)
    if use_database:
        v.OpenDatabase(args.database + ' database', 0)
    else:
        v.OpenDatabase(args.database)

    it = 0
    if not args.i0:
        args.i0 = 0
    if not args.i1:
        args.i1 = v.TimeSliderGetNStates()

    for i in range(args.i0, args.i1):
        v.SetTimeSliderState(i)
        v.AddPlot("Pseudocolor", args.variable)
        p = v.PseudocolorAttributes()
        p.centering = p.Nodal
        #Setting the min and max values of the plot
        if not np.isnan(args.pltLims[0]):
            p.minFlag = 1
            p.min = args.pltLims[0]
        if not np.isnan(args.pltLims[1]):
            p.maxFlag = 1
            p.max = args.pltLims[1]
        if (p.max > 0.0) and (p.min > 0.0) and (args.scale == "Log"):
            p.scaling = p.Log
        else:
            p.scaling = p.Linear
        v.SetPlotOptions(p)
        #Making a reflection of the plot
        v.AddOperator("Reflect", 0)
        ReflectAtts = v.ReflectAttributes()
        ReflectAtts.octant = ReflectAtts.PXPYPZ  # PXPYPZ, NXPYPZ, PXNYPZ, NXNYPZ, PXPYNZ, NXPYNZ, PXNYNZ, NXNYNZ
        ReflectAtts.useXBoundary = 1
        ReflectAtts.specifiedX = 0
        ReflectAtts.useYBoundary = 1
        ReflectAtts.specifiedY = 0
        ReflectAtts.useZBoundary = 1
        ReflectAtts.specifiedZ = 0
        ReflectAtts.reflections = (1, 1, 0, 0, 0, 0, 0, 0)
        ReflectAtts.planePoint = (0, 0, 0)
        ReflectAtts.planeNormal = (0, 0, 0)
        ReflectAtts.reflectType = ReflectAtts.Axis  # Plane, Axis
        v.SetOperatorOptions(ReflectAtts, -1, 0)

        #Settings for the atributes
        atts = v.AnnotationAttributes()
        atts.userInfoFlag = 0
        atts.databaseInfoFlag = 0
        atts.axes2D.visible = 0
        atts.legendInfoFlag = 0
        v.SetAnnotationAttributes(atts)
        #Setting the view attributes
        vat = v.View2DAttributes()
        vat.windowCoords = (-0.02,0.02,0,0.16)
        vat.viewportCoords = (0, 1, 0, 1)
        vat.fullFrameActivationMode = vat.On
        v.SetView2D(vat)
        v.DrawPlots()
        #Saving the figure with settings
        sfatts = v.SaveWindowAttributes()
        sfatts.family = 0
        sfatts.format = sfatts.PNG
        sfatts.fileName = args.database.split("*")[0]+str(i)+"_"+args.variable
        sfatts.outputToCurrentDirectory = 0
        sfatts.outputDirectory = args.outputDir
        #sfatts.screenCapture = 1
        sfatts.resConstraint = sfatts.NoConstraint
        sfatts.width = 1024
        sfatts.height = 3215
        sfatts.quality = 1000
        v.SetSaveWindowAttributes(sfatts)
        print("Saving figure", i)
        v.SaveWindow()
                

    sys.exit()
