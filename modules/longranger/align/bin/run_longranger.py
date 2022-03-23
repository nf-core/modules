#!/usr/bin/env python

from argparse import ArgumentParser
import subprocess
import os.path

def run(sample, fastqs, reference, overrides):
    command = ['longranger', 'align', '--id='+sample, '--fastqs='+fastqs,\
        '--sample='+sample, '--reference='+reference,\
        '--jobmode=lsf', '--localcores=16', '--localmem=60', '--maxjobs=1000',\
        '--jobinterval=100', '--disable-ui', '--nopreflight',\
        '--override='+overrides]
    subprocess.Popen(command) 
    output = sample+'/outs/possorted_bam.bam'
    index = sample+'/outs/summary.csv'
    while (not os.path.isfile(output) or not os.path.isfile(index)) :
        pass 

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('sample')
    parser.add_argument('fastqs')
    parser.add_argument('reference')
    parser.add_argument('overrides')
    args = parser.parse_args()    
    run(args.sample, args.fastqs, args.reference, args.overrides)
