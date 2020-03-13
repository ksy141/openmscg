#!/usr/bin/env python3

import numpy as np
from cg import *

TIMER.reset()

import argparse
parser = argparse.ArgumentParser(description='Run MSCG Range-finder, RDF Calculations, or Inversed-boltzmann method.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("topology",   help="topology file",   type=str)
parser.add_argument("trajectory", help="trajectory file", type=str)

parser.add_argument("--skip",   metavar='N', type=int, default=0, help="skip first N frames")
parser.add_argument("--every",  metavar='N', type=int, default=1, help="read the first frame in every N frames")
parser.add_argument("--frames", metavar='N', type=int, default=0, help="maximum frames to be processed, 0 for all")

parser.add_argument("--top",    metavar='', type=str, default="cgtop", help="type of topology file")
parser.add_argument("--names",  metavar='', type=str, help="comma separated atom names (needed when using LAMMPS data file for topology)")
parser.add_argument("--cut",    metavar='', type=float, default=10.0, help="cut-off for pair interactions")

parser.add_argument("--pair", metavar='', type=str, help="add a pair RDF, in format of [type1,type2,min,res]", action='append')
parser.add_argument("--bond", metavar='', type=str, help="add a bond distribution, in format of [type1,type2,min,max,res]", action='append')
parser.add_argument("--angle", metavar='', type=str, help="add an angle distribution, in format of [type1,type2,type3,min,max,res]", action='append')

args = parser.parse_args()

# Build up topology

top = build_topology(args.top, args.topology, {
    'names': [] if args.names is None else args.names.split(',')
})

top.build_special(True, True, True)

# Load trajectory

trj = Trajectory(args.trajectory)

if top.natoms != trj.natoms:
    raise Exception("Inconsistent number of atoms between topology (%d) and trajectory (%d)." % (top.natoms, trj.natoms))

# Build up lists

plist = PairList(top, trj)
plist.init(cut=args.cut)
plist.setup_bins()
blist = BondList(top, trj)

# Build up tables

if args.pair is not None:
    for pair in args.pair:
        w = pair.split(',')
        TablePairBSpline(plist, top.get_pair_type(w[0], w[1]), xmin=float(w[2]), resolution=float(w[3])).setup_cache()

if args.bond is not None:
    for bond in args.bond:
        w = bond.split(',')
        TableBondBSpline(blist, top.get_bond_type(w[0], w[1]), xmin=float(w[2]), xmax=float(w[3]), resolution=float(w[4])).setup_cache()

if args.angle is not None:
    for angle in args.angle:
        w = angle.split(',')
        TableAngleBSpline(blist, top.get_angle_type(w[0], w[1], w[2]), xmin=float(w[3]), xmax=float(w[4]), resolution=float(w[5])).setup_cache()
    
# Process trajectory

matrix = Matrix()
matrix.add_tables(tables.all)
matrix.setup(top.natoms)

for i in range(args.skip):
    trj.read_frame()

TIMER.click('init')
nread = 0

while trj.read_frame():
    
    TIMER.click('io')
    TIMER.click('matrix', matrix.reset())
    TIMER.click('io', traj.read_frame())
    TIMER.click('pair', plist.build())
    TIMER.click('bond', blist.build())
    TIMER.click('table', tables.compute_all())
    TIMER.click('matrix', matrix.multiplyadd(trj))

    # goto the next frame
    
    for i in range(args.every -1):
        trj.read_frame()
    
    nread += 1
    if args.frames>0 and nread >= args.frames:
        break


matrix.solve()
TIMER.click('solver')

TIMER.report()

# Output


