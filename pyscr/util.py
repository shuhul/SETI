import os
import numpy as np
import math


def find(name):
    for root, dirs, files in os.walk('../'):
        if name in files:
            return os.path.join(root, name)

def load(name, skiprow):
    return np.loadtxt(find(name), skiprows=skiprow)

def log10(vals, default=0):
    return np.log10(vals, out=np.zeros_like(vals), where=(vals>0))

def findAndReplaceLine(file, key, replace):
    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i].find(key) != -1:
                lines[i] = replace
    with open(file, 'w') as f:
        f.writelines(lines)

def yrToSec(yr):
    return 3.154e7*float(yr)

def toSciNot(num):
    return '{:e}'.format(num)




