#!/usr/bin/env python2.7

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import visit as v
import sys


def get_args():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Save lineout data using Visit in curve format.
        The output file contains two columns: the path length and the
        variable along the path.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s lineout.py database var
        -r0 x y z -r1 x y z''')
    pr.add_argument('database', type=str,
                    help='Database name (e.g. sim.silo or "sim_*.silo")')
    pr.add_argument('varname', type=str,
                    help='Extract this (scalar) variable')
    pr.add_argument('-r0', type=float, nargs=3, required=True,
                    metavar=('x0', 'y0', 'z0'),
                    help='Start location of line (z=0 in 2D)')
    pr.add_argument('-r1', type=float, nargs=3, required=True,
                    metavar=('x1', 'y1', 'z1'),
                    help='End location of line (z=0 in 2D)')
    pr.add_argument('-outdir', type=str, default=os.getcwd(),
                    help='Output directory')
    pr.add_argument('-fname', type=str, default='line',
                    help='Output filename (without extension)')
    return pr.parse_args()


if __name__ == '__main__':
    args = get_args()

    # Open either a database or a single file
    if '*' in args.database:
        print("database found")
        v.OpenDatabase(args.database + ' database', 0)
    else:
        v.OpenDatabase(args.database)

    v.AddPlot("Curve", "operators/Lineout/" + args.varname, 1, 1)

    LineoutAtts = v.LineoutAttributes()
    LineoutAtts.point1 = tuple(args.r0)
    LineoutAtts.point2 = tuple(args.r1)
    LineoutAtts.interactive = 0
    LineoutAtts.ignoreGlobal = 0
    v.SetOperatorOptions(LineoutAtts, 1)

    v.DrawPlots()

    # Set output options
    s = v.SaveWindowAttributes()
    s.format = s.CURVE
    s.fileName = args.fname + '_' + args.varname + '_'
    s.outputDirectory = args.outdir
    s.outputToCurrentDirectory = 0
    v.SetSaveWindowAttributes(s)

    imax = v.TimeSliderGetNStates()
    imin = 0

    for i in range(imin, imax):
        v.TimeSliderSetState(i)
        v.SaveWindow()

    sys.exit()
