import time
import matplotlib.pyplot as plt
import numpy as np

def Performance(functions, dataset, figsize=plt.figaspect(1), 
                N=3, logx=True, logy=True, out=None,
                fontsize=14.0, fontfamily='Arial', mathfont='stix'):
    
    plt.rcParams['font.size']        = fontsize
    plt.rcParams['font.family']      = fontfamily
    plt.rcParams['mathtext.fontset'] = mathfont

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlabel('size of data')
    ax.set_ylabel('time')

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    if not isinstance(functions, list):
        functions = [functions]
    if not isinstance(dataset, list):
        dataset   = [dataset]
    

    for function in functions:
        xdata     = []
        ydata_avg = []
        ydata_std = []
        for data in dataset:
            xdata.append(len(data))
            
            y = []
            for i in range(N):
                t1 = time.time()
                function(data) 
                t2 = time.time()
                y.append(t2 - t1)

            ydata_avg.append(np.average(y))
            ydata_std.append(np.std(y))

        ax.errorbar(xdata, ydata_avg, yerr=ydata_std)
    
    fig.tight_layout()
    if out:
        fig.savefig(out)
    else:
        plt.show()
            

