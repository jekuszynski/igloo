#################################################################
#                                                               #
# Credit to https://www.jappoker.com/blog/2019/band-diagram-QE/ #
# for the original code this is based on.                       #
#                                                               #
#################################################################
import numpy as np
import matplotlib.pyplot as plt
import os
import re

def dist(a,b):
    # calculate the distance between vector a and b
    d = ((a[0]-b[0]) ** 2 + (a[1]-b[1]) ** 2 + (a[2]-b[2]) ** 2)**0.5
    return d

def read_fermi(file_name):
    # Read the fermi energy in scf.out
    fermi = 0
    with open(file_name, "r") as f:
        lines = f.readlines()
    for line in lines:
        if "the Fermi energy" in line:
            fermi = float(line.split()[4])
    
    return fermi

def read_bnd(file_name):
    # Read the bands in Band.dat
    coord_regex = r"^\s+(.\d+.\d+)\s+(.\d+.\d+)\s+(.\d+.\d+)$"
    x_coord = []
    x = []
    bands = dict()

    with open(file_name, "r") as f:
        lines = f.readlines()

    for i in range(len(lines)):
        line = lines[i]
        match = re.match(coord_regex,line)
        if match:
            x_coord.append([float(match.group(1)), float(match.group(2)), float(match.group(3)) ])
            bandddd = lines[i+1] + lines[i+2]
            bandddd = bandddd.split()

            for j in range(len(bandddd)):
                if j not in bands.keys():
                    bands[j] = []
                bands[j].append(float(bandddd[j]))
    for i in range(len(x_coord)) :
        if i == 0:
            x.append(0)
        else:
            x.append( x[-1] + dist(x_coord[i], x_coord[i-1]))
    return bands,x

def plot(bands, x, fermi):
    xaxis = [min(x),max(x)]
    for i in bands.values(): 
        plt.plot(x, i, color=color_dic[pressure], lw=0.2)
    plt.plot(xaxis, [fermi, fermi], color="#66ccff",ls="solid", alpha = 0.5,lw = 1.2)
    plt.xlim(xaxis)