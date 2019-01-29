#!/usr/bin/python3

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-r", type=float, default=None, help="radius")
parser.add_argument("-g", type=float, default=None, help="gamma")
parser.add_argument("-m", type=int,   default=None, help="mobile particles")
parser.add_argument("--sigma", type=float, default=1.0, help="lennard jones sigma")
args = parser.parse_args()

r_opt = np.power(2.0, 1.0/6.0) * args.sigma

if args.r != None and args.g != None and args.m != None:
    raise ValueError("please choose 1 input parameters, not all 3")

if args.r != None and args.g != None and args.m == None:
    raise ValueError("please choose 1 input parameters, not 2")
elif args.r != None and args.g == None and args.m != None:
    raise ValueError("please choose 1 input parameters, not 2")
elif args.r == None and args.g != None and args.m != None:
    raise ValueError("please choose 1 input parameters, not 2")

if args.r != None:
    args.g = np.arcsin(r_opt/(2.0*args.r))
    args.m = 16.0*args.r*args.r / (1.1027*r_opt*r_opt)
elif args.g != None:
    args.g *= np.pi/180
    args.r = r_opt / (2.0*np.sin(args.g))
    args.m = 4.0 / (1.1027*np.sin(args.g)*np.sin(args.g))
elif args.m != None:
    args.r = r_opt*np.sqrt(1.1027*args.m) / 4.0
    args.g = np.arcsin(2.0 / (np.sqrt(1.1027*args.m)))

print(f"radius: {args.r:>.4f}   gamma: {args.g*180/np.pi:>.4f}   particles: {args.m:>.1f}")