import os
import sys
import argparse


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser(description='QC analysis')
    parser.add_argument('--inputfiles', required=True, nargs='+')
    parser.add_argument('--outputdir', required=True)
    parser.add_argument('--ignorefirst', default=600, type=int)
    parser.add_argument('--halfwindow', default=10, type=int)
    parser.add_argument('--peakthreshold', default=3, type=float)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # make output dir
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    # loop over input files
    cmds = []
    for inputfile in args.inputfiles:

        # define output file
        base = os.path.basename(inputfile).replace('.csv','.png')
        outputfile = os.path.join(args.outputdir,base)

        # make command
        cmd = 'python3 qcanalysis.py'
        cmd += ' --inputfile {}'.format(inputfile)
        cmd += ' --outputfile {}'.format(outputfile)
        cmd += ' --ignorefirst {}'.format(args.ignorefirst)
        cmd += ' --halfwindow {}'.format(args.halfwindow)
        cmd += ' --peakthreshold {}'.format(args.peakthreshold)
        cmds.append(cmd)

    # run commands
    for cmd in cmds: os.system(cmd)
