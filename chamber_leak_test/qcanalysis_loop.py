import os
import sys
import argparse


if __name__=='__main__':

    # read command line arguments
    parser = argparse.ArgumentParser(description='QC analysis')
    parser.add_argument('--inputfiles', required=True, nargs='+')
    parser.add_argument('--outputdir', required=True)
    parser.add_argument('--ignorefirst', default=600, type=int)
    args = parser.parse_args()

    # print arguments
    print('Running with following configuration:')
    for arg in vars(args):
        print('  - {}: {}'.format(arg,getattr(args,arg)))

    # make output dir
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    # group input files in bottom and top gap files
    bases = []
    for f in args.inputfiles:
        base = os.path.basename(f).rsplit('-',1)[0]
        bases.append(base)
    bases = sorted(list(set(bases)))
    cfiles = []
    for base in bases:
        thiscfiles = [base, None, None]
        for f in args.inputfiles:
            if os.path.basename(f).startswith(base) and 'BOTTOM' in f:
                thiscfiles[1] = f
            if os.path.basename(f).startswith(base) and 'TOP' in f:
                thiscfiles[2] = f
        if None in thiscfiles:
            print('WARNING: no matching files found for base {}'.format(base))
        else:
            cfiles.append(thiscfiles)

    # loop over input files
    cmds = []
    for thiscfiles in cfiles:

        # define output file
        outputfile = os.path.join(args.outputdir,thiscfiles[0]+'.png')

        # make command
        cmd = 'python3 qcanalysis.py'
        cmd += ' --bottomgapfile {}'.format(thiscfiles[1])
        cmd += ' --topgapfile {}'.format(thiscfiles[2])
        cmd += ' --outputfile {}'.format(outputfile)
        cmd += ' --bottom_ignorefirst {}'.format(args.ignorefirst)
        cmd += ' --top_ignorefirst {}'.format(args.ignorefirst)
        cmds.append(cmd)

    # run commands
    for cmd in cmds:
        print(cmd)
        os.system(cmd)
