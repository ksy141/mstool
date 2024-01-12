import numpy as np
import argparse
import matplotlib.pyplot as plt
from   matplotlib.animation import FuncAnimation
from   matplotlib.lines     import Line2D

description = "This script gets timeseries data and makes an animation plot. Always needs to provide xlim ylim."

def parse_args():
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Positional arguments
    parser.add_argument("in_csv",  
                        type  = str, 
                        help  = 'data file (e.g., csv)')

    parser.add_argument("out_video", 
                        type  = str, 
                        help  = 'video file (e.g., mp4)')

    # Optional arguments
    parser.add_argument("--fps",        
                        type    = float,
                        default = 25, 
                        help    = 'frame per second in the output')

    parser.add_argument("--fontsize",   
                        type    = float,
                        default = 14.0)

    parser.add_argument("--fontfamily", 
                        default = 'Arial')

    parser.add_argument("--mathtext",   default='stix')
    parser.add_argument("--xlabel",     default=None)
    parser.add_argument("--ylabel",     default=None)
    parser.add_argument("--sep",        default=",")
    parser.add_argument("--comment",    default="#")
    parser.add_argument("--skiprows",   default=0,   type=int)
    parser.add_argument("--dpi",        default=300, type=int)
    parser.add_argument("--xtimes",     default=1.0, type=float)
    parser.add_argument("--ytimes",     default=1.0, type=float)

    parser.add_argument("--figsize", 
                        nargs   = '*', 
                        default = [5,5],
                        type    = float)

    parser.add_argument("--xlim", 
                        nargs   = '*', 
                        default = [None],
                        type    = float)

    parser.add_argument("--ylim", 
                        nargs   = '*', 
                        default = [None],
                        type    = float)

    parser.add_argument("--xticks", 
                        nargs   = '*', 
                        default = [None],
                        type    = float)

    parser.add_argument("--yticks", 
                        nargs   = '*', 
                        default = [None],
                        type    = float)

    parser.add_argument("--marker",
                        nargs   = '*',
                        default = [None],
                        help    = 'format to depict data (e.g., C0 ro)')

    parser.add_argument("--linewidth", "--lw",
                        nargs   = '*',
                        default = [1.5])

    parser.add_argument("--markersize", "--ms",
                        nargs   = '*',
                        default = [None])

    parser.add_argument("--xindex", 
                        default = 0,
                        type    = int, 
                        help    = "the column number of the x (0-indexed) of the input file")

    parser.add_argument("--yindex",
                        nargs   = '*',
                        default = [1],
                        type    = int, 
                        help    = "the column number of the y (0-indexed) of the input file")

    parser.add_argument("--bar",
                        action  = 'store_true',
                        dest    = 'bar',
                        help    = 'plot data first and move a bar')
    parser.add_argument("--no-bar",
                        action  = 'store_false',
                        dest    = 'bar',
                        help    = 'do not use a bar (default)')
    parser.set_defaults(bar=False)

    parser.add_argument("--barlinecolor",
                        default = 'black',
                        help    = 'line color of the bar')

    parser.add_argument("--barlinestyle",
                        default = '-',
                        help    = 'line style of the bar')


    parser.add_argument("--barlinewidth",
                        default = 1.5,
                        type    = float,
                        help    = 'linewidth of the moving bar')

    parser.add_argument("--barmarkersize",
                         default = None,
                         type    = float,
                         help    = "size of marker at intersection of bar and data")

    parser.add_argument("--barmarker",
                         default = 'o',
                         help    = "marker at intersection of bar and data")

    parser.add_argument("--startindex",
                         default = None,
                         type    = int,
                         help    = "starting index (1-index; to be consistent with ChimeraX) of bar")

    parser.add_argument("--endindex",
                         default = None,
                         type    = int,
                         help    = "ending index (1-index; to be consistent with ChimeraX; including this index) of bar")

    parser.add_argument("--black", 
                        action  = 'store_true', 
                        dest    = 'black', 
                        help    = 'black background')
    parser.add_argument('--no-black',
                        action  = 'store_false',
                        dest    = 'black',
                        help    = 'white background (default)')
    parser.set_defaults(black=False)

    return parser.parse_args()



def update_lines(num, data, lines):
    for i, line in enumerate(lines):
        line.set_data(data[:num,0], data[:num,i+1])
    return lines

def update_bar(num, data, lines, ylim):
    for i, line in enumerate(lines):
        # i = 0 -> bar
        x = data[num, 0]
        y = data[num, i]

        if isinstance(line, Line2D):
            # bar
            #line.set_xdata([x])
            line.set_xdata([x, x])
            line.set_ydata(ylim)

        else:
            # scatter
            line.set_offsets([x, y])

    return lines


def setting_params(value, yindex):
    if len(value) == len(yindex):
        vals = value
    elif len(value) == 1:
        vals = value * len(yindex)
    else:
        raise ValueError('Either len(value) should be len(args.yindex) or 1')
    return vals
     

def main():
    ### PARSING
    args = parse_args()

    ### SETTING MATPLOTLIB
    plt.rcParams['font.size']        = args.fontsize
    plt.rcParams['font.family']      = args.fontfamily
    plt.rcParams['mathtext.fontset'] = args.mathtext

    ### BLACK?
    if args.black:
        plt.rcParams.update({
            "lines.color": "white",
            "patch.edgecolor": "white",
            "text.color": "white",
            "axes.facecolor": "black",
            "axes.edgecolor": "lightgray",
            "axes.labelcolor": "white",
            "xtick.color": "white",
            "ytick.color": "white",
            "grid.color": "lightgray",
            "figure.facecolor": "black",
            "figure.edgecolor": "black",
            "savefig.facecolor": "black",
            "savefig.edgecolor": "black"})


    ### READ DATA
    rawdata = []
    with open(args.in_csv) as fp:
        for i, line in enumerate(fp):
            if i < args.skiprows: continue
            if line.startswith(args.comment): continue
            sl = line.strip().split(args.sep)
            rawdata.append([float(c) for c in sl])
    rawdata = np.array(rawdata)
    data    = rawdata[:, [args.xindex, *args.yindex]]

    data[:,0]  *= args.xtimes
    data[:,1:] *= args.ytimes
    
    ### MAKE FIGURE
    fig, ax = plt.subplots(figsize=args.figsize, dpi=300)
    if args.xlabel: ax.set_xlabel(args.xlabel)
    if args.ylabel: ax.set_ylabel(args.ylabel)
    if args.xticks[0] != None: ax.set_xticks(args.xticks)
    if args.yticks[0] != None: ax.set_yticks(args.yticks)

    if args.xlim[0] == None: 
        ax.set_xlim(data[:,0].min(), data[:,0].max())
    else:
        ax.set_xlim(args.xlim)

    if args.ylim[0] == None:
        ax.set_ylim([data[:,1:].min(), data[:,1:].max()])
    else:
        ax.set_ylim(args.ylim)
    
    ### LINE COLORS
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors     = prop_cycle.by_key()['color']
    if args.marker[0] == None:
        colors  = colors * (len(args.yindex) // 10 + 1)
        markers = colors[:len(args.yindex)]
    else:
        markers = setting_params(args.marker, args.yindex)

    linewidths  = setting_params(args.linewidth,  args.yindex)
    markersizes = setting_params(args.markersize, args.yindex)
    
    if args.bar:
        # plot data
        colors = []
        for ii, (marker, linewidth, markersize) in enumerate(zip(markers, linewidths, markersizes), 1):
            tmpax = ax.plot(data[:,0], data[:,ii], marker, lw=linewidth, ms=markersize)
            colors.append(tmpax[0].get_color())
        ylim = ax.get_ylim()

        # make empty bar line
        #lines = [ax.axvline([], c=args.barlinecolor, ls=args.barlinestyle, lw=args.barlinewidth)]
        lines = [ax.plot([], [], lw=args.barlinewidth, ls=args.barlinestyle, c=args.barlinecolor)[0]]
        
        # make empty scatter data
        for color in colors:
            lines.append(ax.scatter([], [], s=args.barmarkersize, color=color, marker=args.barmarker))

        # make animation
        if isinstance(args.startindex, int):
            startindex = args.startindex - 1
        else:
            startindex = None

        if isinstance(args.endindex, int):
            # -1 for 1-index; +1 for including the last value
            endindex   = args.endindex - 1 + 1
        else:
            endindex   = None

        xrange = [i for i in range(len(data))][startindex:endindex]
        ani = FuncAnimation(fig, update_bar, xrange, fargs=(data, lines, ylim), blit=True)

    else:
        # make empty
        lines = []
        for marker, linewidth, markersize in zip(markers, linewidths, markersizes):
            lines.append(ax.plot([], [], marker, lw=linewidth, ms=markersize)[0])

        # make animation
        ani = FuncAnimation(fig, update_lines, len(data), fargs=(data, lines), blit=True)
    
    fig.tight_layout()
    ani.save(args.out_video, fps=args.fps)


if __name__ == '__main__':
    main()

