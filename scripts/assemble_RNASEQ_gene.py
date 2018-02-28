#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import glob, os
import argparse
import textwrap
from os import listdir
from os.path import isfile,join
from sys import argv


all_rnaseq_files = sorted([f for f in listdir(argv[1]) if (f.endswith('genes.results') and isfile(join(argv[1],f)))])
## print(all_rnaseq_files)

samplename = []
tpm_matrix = []
samplename.append("Gene")
for f in all_rnaseq_files:
    samplename.append(f.split('.genes.results')[0])
    genename = []
    tpm = []
    rnaseq_file = open(argv[1] + f,'r')
    rnaseq_file_lines = rnaseq_file.readlines()
    header = rnaseq_file_lines[0]
    tpm_index = header.strip().split('\t').index("TPM")
    for line in rnaseq_file_lines[1:]:
        genename.append(line.rstrip("\r\n").split('\t')[0])
        tpm.append(line.rstrip("\r\n").split('\t')[tpm_index])
    tpm_matrix.append(tpm)

#tpm_matrix = zip(*tpm_matrix)
output_matrix = zip(genename,*tpm_matrix)
#print(output_matrix[0:10])



with open(argv[2],'w') as f:
   writer = csv.writer(f,delimiter=',')
   writer.writerow(samplename)
   writer.writerows(output_matrix)
