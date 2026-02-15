import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


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


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser(description='QC analysis')
    parser.add_argument('-b', '--bottomgapfile', required=True)
    parser.add_argument('-t', '--topgapfile', required=True)
    parser.add_argument('--bottom_ignorefirst', default=600, type=int)
    parser.add_argument('--bottom_ignorelast', default=0, type=int)
    parser.add_argument('--top_ignorefirst', default=600, type=int)
    parser.add_argument('--top_ignorelast', default=0, type=int)
    parser.add_argument('-o', '--outputfile', default=None)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # get the data
    bottomdf = pd.read_csv(args.bottomgapfile, sep=',', header=0)
    bottomdata = bottomdf.values
    topdf = pd.read_csv(args.topgapfile, sep=',', header=0)
    topdata = topdf.values

    # define inline plotting function
    def do_fit(data, ignorefirst=0, ignorelast=0, degree=1):
        linwindow = data[ignorefirst:-1-ignorelast,0]
        lindata = nan_to_num(data[ignorefirst:-1-ignorelast,1])
        try: params = np.polyfit(linwindow, lindata, degree)[::-1]
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
        pdrop = -(params[1])*600
        pdropstr = 'Pressure drop: {:.2f} mbar / 10 mins'.format(pdrop)
        pdropcolor = 'darkgreen' if pdrop<0.4 else 'red'
        return {'xax': linwindow, 'fit': linfit,
                'pdropstr': pdropstr, 'pdropcolor': pdropcolor}

    # do the fits
    bottomfitresults = do_fit(bottomdata, 
            ignorefirst=args.bottom_ignorefirst,
            ignorelast=args.bottom_ignorelast,
            degree=1)
    topfitresults = do_fit(topdata,
            ignorefirst=args.top_ignorefirst,
            ignorelast=args.top_ignorelast,
            degree=1)

    # make plot
    fig,axs = plt.subplots(figsize=(6,8), nrows=2, sharex=False)
    plt.subplots_adjust(hspace=0.2)
    for ax in axs:
        #ax.grid(visible=True)
        ax.tick_params(direction='in', which='both')

    # plot of top gap
    axs[0].plot(topdata[:,0], topdata[:,1], color='mediumblue', label='Raw data')
    axs[0].plot(topfitresults['xax'], topfitresults['fit'],
            color='red', label = 'Linear fit')
    (ymin,ymax) = axs[0].get_ylim()
    axs[0].set_ylim(ymin-(ymax-ymin)*0.1, ymax+(ymax-ymin)*0.2)
    (newymin,newymax) = axs[0].get_ylim()
    axs[0].vlines([args.top_ignorefirst, len(topdata)-args.top_ignorelast],
            ymin, ymin + 0.8*(ymax-ymin),
            color='grey', linestyles='--')
    axs[0].legend()
    pdroptxt = axs[0].text(0.95, 0.75, topfitresults['pdropstr'],
            horizontalalignment='right',
            transform=axs[0].transAxes,
            color=topfitresults['pdropcolor'])
    axs[0].set_ylabel('Pressure (mbar)')
    axs[0].set_xlabel('Time (s)')
    axs[0].text(0.05, 0.05, 'TOP GAP', transform=axs[0].transAxes)

    # plot of bottom gap
    axs[1].plot(bottomdata[:,0], bottomdata[:,1], color='mediumblue', label='Raw data')
    axs[1].plot(bottomfitresults['xax'], bottomfitresults['fit'],
            color='red', label = 'Linear fit')
    (ymin,ymax) = axs[1].get_ylim()
    axs[1].set_ylim(ymin-(ymax-ymin)*0.1, ymax+(ymax-ymin)*0.2)
    (newymin,newymax) = axs[1].get_ylim()
    axs[1].vlines([args.bottom_ignorefirst, len(bottomdata)-args.bottom_ignorelast],
            ymin, ymin + 0.8*(ymax-ymin), 
            color='grey', linestyles='--')
    axs[1].legend()
    pdroptxt = axs[1].text(0.95, 0.75, bottomfitresults['pdropstr'],
            horizontalalignment='right',
            transform=axs[1].transAxes, 
            color=bottomfitresults['pdropcolor'])
    axs[1].set_ylabel('Pressure (mbar)')
    axs[1].set_xlabel('Time (s)')
    axs[1].text(0.05, 0.05, 'BOTTTOM GAP', transform=axs[1].transAxes)

    # make a title
    title = args.bottomgapfile.rsplit('-',1)[0]
    title = os.path.basename(title)
    if os.path.basename(args.topgapfile.rsplit('-',1)[0])!=title:
        print('WARNING: titles derived from bottom gap file and top gap file'
                +' do not seem to correspond; are they for the same chamber?')
    axs[0].set_title(title, loc='left')

    # show or save the plot
    if args.outputfile is None: plt.show()
    else: fig.savefig(args.outputfile)
