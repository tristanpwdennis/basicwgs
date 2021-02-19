#basic system wrapper for processing and qc of wgs data
#in the absence of a system for using nextflow
#tristan dennis feb 2021


import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("ref", help="specify the location of your reference assembly", type=str)
args = parser.parse_args()
print(args.ref)
