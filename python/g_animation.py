"""
Module to plot Heavy output

Useful data are generated on a grid of stations ("grid" output, i.e., -g flag)

Generates two plots:
1) map-view animation of gravity change, saved as a .gif file
2) time-series plot of gravity change at all locations (shown after closing 
   first plot window)

Requires Heavy .lst and .out files. Both must have the same basename (e.g., 
model1.lst and model1.out).

Example:
$> python g_animation.py model1

positional arguments:
  prefix      Heavy prefix for .lst and lout files

optional arguments:
  -h, --help            show help message and exit
  -v, --version         show program's version number and exit
  -l LAYER, --layer LAYER
                        Plot a specific layer
"""

import os, sys
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import matplotlib.animation as animation
import argparse

HEADER_STRING = 'HEAVY-GENERATED GRAVITY OBSERVATION LOCATIONS'

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

def plot_ts(g_data):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # for row in g_data:
    ax.plot(g_data)
    ax.set_ylabel(r'Gravity change, uGal')
    plt.show()

def init_plot(args):
    lbl, x, y, z = load_coords(args.prefix[0])
    g_lbl, ts, data = load_data(args.layer, args.prefix[0])
    unique_ts = list(set(ts))
    unique_ts.sort(key=float)
    g_data = []
    for i, this_ts in enumerate(unique_ts):
        g_data.append(get_ts_data(this_ts, ts, data))
    g_data = np.array(g_data)
    animate_heavy(unique_ts, g_data, x, y, z)
    plot_ts(g_data)

def decimate_data(x, y, z, data):
    # Remove points at 0 elevation (outside model bounds)
    xp, yp, dp, idxs  = [], [], [], []
    for idx in range(len(z)):
        if float(z[idx]) > 0.001:
            xp.append(float(x[idx]))
            yp.append(float(y[idx]))
            # dp.append(float(data[idx]))
            idxs.append(idx)
    dp = [float(d) for d in data]
    return xp, yp, dp, idxs

def animate_heavy(unique_ts, g_data, x, y, z):
    xp, yp, gd, idxs = decimate_data(x, y, z, g_data[0])
    scat = plt.scatter(xp, yp, c=gd, s=20)
    scat.set_clim(get_clim(g_data))  # [g_data[0].min(),10])
    plt.gca().set_aspect('equal')
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    a = plt.colorbar()
    a.set_label('Gravity change, in uGal')
    ani = animation.FuncAnimation(fig2, update_plot, frames=range(len(unique_ts)),
                                  fargs=(g_data, idxs, scat, unique_ts), interval=400)
    ani.save('heavy_output.gif')
    plt.show()
    
def get_clim(data):
    mn = min(min(data, key=min))    
    mx = max(max(data, key=max))    
    return [mn, mx]
    
def update_plot(i, data, idxs, scat, ts):
    dp = []
    d = data[i]
    scat.set_array(np.array(d))
    #scat.set_clim([d.min(), d.max()])
    ax2.set_title("{:.0f} years".format(float(ts[i])/365))
    plt.draw()
    labels = [item.get_text() for item in ax2.get_xticklabels()]

def get_ts_data(this_ts, ts, data):
    out = []
    for i, row in enumerate(ts):
        if row == this_ts:
            out.append(data[i])
    return out

def load_data(layer, fn):
    lbl, ts, g = [], [], []
    with open(fn + '.out', 'r') as fid:
        for line in fid:
            line_elems = line.split()
            lbl.append(line_elems[0])
            ts.append(line_elems[3])
            try:
                if layer:
                    g.append(float(line_elems[3+int(layer)]))
                else:
                    g.append(float(line_elems[-1]))
            except ValueError as e:
                print("Error reading file. This routine doesn't work unless the -g option was used with Heavy.")
                sys.exit(1)
    return (lbl, ts, g)

def load_coords(fn):
    g_lbl, g_x, g_y, g_z = [], [], [], []
    with open(fn + '.lst', 'r') as fid:
        for line in fid:
            if HEADER_STRING in line:
                break
        _ = fid.readline()
        _ = fid.readline()
        line = fid.readline().split()
        while line:
            g_lbl.append(line[0])
            g_z.append(line[1])
            g_x.append(line[2])
            g_y.append(line[3])
            line = fid.readline().split()
    return (g_lbl, g_x, g_y, g_z)

def init_argparse():
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Plot HEAVY results"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = "{} version 1.0.0".format(parser.prog)
    )
    parser.add_argument(
        "-l", "--layer", action="store",
        help="Plot a specific layer")
    parser.add_argument('prefix', help='Heavy prefix for .lst and lout files', nargs='*')
    return parser

if __name__ == '__main__':
    fn = 'tran'
    parser = init_argparse()
    args = parser.parse_args()
    init_plot(args)
