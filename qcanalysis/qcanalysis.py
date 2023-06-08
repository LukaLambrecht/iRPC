import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


def smooth(series, halfwindow=1):
    ksize = 2*halfwindow+1
    kernel = np.ones(ksize)/ksize
    smooth = np.convolve(series, kernel, mode='same')
    for i in range(halfwindow):
        smooth[i] = series[halfwindow]
        smooth[-1-i] = series[-1-halfwindow]
    return smooth

def nan_to_num(series):
    naninds = np.argwhere(np.isnan(series))
    startidx = 1
    endidx = 0
    res = np.copy(series)
    while startidx<len(series)-1:
        if not np.isnan(series[startidx]):
            startidx += 1
            continue
        endidx = startidx + 1
        while endidx<len(series):
            if np.isnan(series[endidx]):
                endidx += 1
                continue
            else: break
        startvalue = series[startidx-1]
        endvalue = series[endidx]
        for idx in range(startidx, endidx):
            intervalue = startvalue + (endvalue-startvalue)/(endidx-startidx)*(idx-startidx)
            res[idx] = intervalue
        startidx = endidx
        continue
    return res

def find_threshold(dev):
    lowthreshold = np.quantile(dev, 0.1)
    highthreshold = np.quantile(dev, 0.9)
    return (lowthreshold+highthreshold)*0.75


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser(description='QC analysis')
    parser.add_argument('--inputfile', required=True)
    parser.add_argument('--ignorefirst', default=600, type=int)
    parser.add_argument('--halfwindow', default=10, type=int)
    parser.add_argument('--peakthreshold', default=3, type=float)
    parser.add_argument('--outputfile', default=None)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # get the data
    df = pd.read_csv(args.inputfile, sep=',', header=0)
    data = df.values

    # find start of fluctuations
    smoothed = smooth(data[:,1], halfwindow=args.halfwindow)
    smoothed = nan_to_num(smoothed)
    dev = np.abs(smoothed-nan_to_num(data[:,1]))
    devsmoothed = smooth(dev, args.halfwindow)
    for i in range(2*args.halfwindow):
        devsmoothed[i] = 0
        devsmoothed[-1-i] = 0
    devthreshold = find_threshold(devsmoothed)
    fluctvalid = True
    try: fluctidx = np.nonzero(devsmoothed>devthreshold)[0][0] - 2*args.halfwindow
    except:
        print('WARNING: could not determine index of fluctuation start, will guess...')
        fluctidx = args.ignorefirst + 600
        fluctidx = min(fluctidx,len(data))
        fluctvalid = False
    if fluctidx < args.ignorefirst:
        print('WARNING: index of fluctuation start is below ingore window, will guess...')
        fluctidx = args.ignorefirst + 600
        fluctidx = min(fluctidx,len(data))
        fluctvalid = False

    # do linear fit in window
    linwindow = data[args.ignorefirst:fluctidx,0]
    lindata = nan_to_num(data[args.ignorefirst:fluctidx,1])
    try: params = np.polyfit(linwindow, lindata, 1)[::-1]
    except:
        print('WARNING: could not perform linear fit.')
        params = [0,0]
    linfit = np.zeros(len(linwindow))
    fitstr = ''
    for degree,param in enumerate(params):
        linfit += param * np.power(linwindow,degree)
        strfrag = '* x$^{}$'.format(degree)
        if degree==0: strfrag = ''
        if degree==1: strfrag = '* x'
        fitstr += ' + {:.1e} {}'.format(param,strfrag)
    fitstr = fitstr.strip(' +')
    fitstr = 'Fit: y = ' + fitstr
    pdrop = abs(params[1])*600
    pdropstr = 'Pressure drop: {:.2f} mbar / 10 mins'.format(pdrop)
    pdropcolor = 'darkgreen' if pdrop<0.4 else 'red'

    # find peaks
    peaks = []
    if fluctvalid:
        fluctdevsmoothed = devsmoothed[fluctidx:]
        peakthreshold = np.quantile(fluctdevsmoothed,0.5)*args.peakthreshold
        peakinds = np.where(fluctdevsmoothed>peakthreshold)[0]
        thispeak = []
        for i,indx in enumerate(peakinds):
            if( i==0 ):
                thispeak.append(indx)
                continue
            if( indx-1 in thispeak ):
                thispeak.append(indx)
            if( i==len(peakinds)-1 or indx-1 not in peakinds ):
                peaks.append(int(np.mean(thispeak)))
                thispeak = [indx]
                continue
        for i,indx in enumerate(peaks):
            peaks[i] = fluctidx + indx

    # make plot
    fig,axs = plt.subplots(figsize=(6,8), nrows=3, sharex=True,
            gridspec_kw={'height_ratios':(0.5,0.5,1)})
    plt.subplots_adjust(hspace=0.05)
    for ax in axs:
        #ax.grid(visible=True)
        ax.tick_params(direction='in', which='both')

    # first plot of raw data and smoothing
    axs[0].plot(data[:,0], data[:,1], color='mediumblue', label='Raw data')
    axs[0].plot(data[:,0], smoothed, color='fuchsia', label='Smoothed')
    axs[0].legend()
    axs[0].set_ylabel('Pressure (mbar)')
    title = args.inputfile.replace('.csv','')
    title = os.path.basename(title)
    axs[0].set_title(title, loc='left')
    
    # second plot of deviation and thresholds
    axs[1].plot(data[:,0], devsmoothed, color='mediumblue', label='Deviation from smoothed')
    axs[1].plot(data[:,0], np.ones(len(dev))*devthreshold, color='darkviolet',
            linestyle='--', label='Threshold for fluctuation')
    if fluctvalid: axs[1].plot(data[:,0], np.ones(len(dev))*peakthreshold, color='red', 
                     linestyle='--', label='Threshold for peak')
    axs[1].legend()
    axs[1].set_ylim(0, axs[1].get_ylim()[1]*1.3)
    axs[1].set_ylabel('Pressure (mbar)')

    # third plot with linear fit and peaks
    axs[2].plot(data[:,0], data[:,1], color='mediumblue', label='Raw data')
    axs[2].plot(linwindow, linfit, color='red', label = 'Linear fit')
    (ymin,ymax) = axs[2].get_ylim()
    axs[2].set_ylim(ymin-(ymax-ymin)*0.1, ymax+(ymax-ymin)*0.2)
    (newymin,newymax) = axs[2].get_ylim()
    axs[2].vlines([args.ignorefirst,fluctidx], ymin, ymin + 0.8*(ymax-ymin), 
            color='grey', linestyles='--')
    axs[2].legend()
    for peak in peaks:
        axs[2].arrow(peak, newymin, 0, (ymin-newymin),
                head_length=(ymin-newymin)*0.3, head_width=25, edgecolor='r', facecolor='r')
    pdroptxt = axs[2].text(0.95, 0.75, pdropstr, horizontalalignment='right',
            transform=axs[2].transAxes, color=pdropcolor)
    axs[2].set_ylabel('Pressure (mbar)')
    axs[2].set_xlabel('Time (s)')

    # show or save the plot
    if args.outputfile is None: plt.show()
    else: fig.savefig(args.outputfile)
