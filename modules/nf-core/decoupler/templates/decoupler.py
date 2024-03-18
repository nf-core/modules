#!/usr/bin/env python

import decoupler as dc
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("mat", type=str,
                    help="Text file containing samples as rows and features as columns")
parser.add_argument("net", type=str,
                    help="network in long format")
parser.add_argument("method", choices=['aucell', 'gsea', 'gsva', 'mdt', 'mlm', 'ora', 'udt', 'ulm', 'viper', 'wmean', 'wsum'])

parser.add_argument("--args", type=dict,
                    help='method-specific arguments')

args = parser.parse_args()

mat = pd.read_csv(args.mat, sep='\t')
net = pd.read_csv(args.net, sep='\t')

act_estimates, pvalues = dc.decouple(
    mat=mat,
    net=net,
    methods=args.method,
    args=args.args
)

dc_output = pd.concat([act_estimates, pvalues], axis=0)
dc_output.to_csv(output_file, sep='\t')

## VERSIONS FILE
with open('versions.yml', 'a') as version_file:
    version_file.write('"${task.process}"\n')
    version_file.write("\tdecoupler: " + dc.__version__ + '\n')
