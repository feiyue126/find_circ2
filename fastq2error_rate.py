#!/usr/bin/env python3

import sys
import math

__version__ = "0.1"
__author__ = "Marcel Schilling"
__credits__ = ["Marvin Jens","Marcel Schilling","Petar Glazar","Nikolaus Rajewsky"]
__status__ = "beta"
__licence__ = "GPL"
__email__ = "marcel.schilling@mdc-berlin.de"

lines_per_fastq_record=4
phred_shift=33
phred_factor=-10
phred_base=10.0

transform_base=math.pow(phred_base,1/phred_factor)
transform_scale=math.pow(transform_base,phred_shift)

lines_since_phred_line=0
nc_count=0
mean=0

for line in sys.stdin:
    lines_since_phred_line=(lines_since_phred_line+1) % lines_per_fastq_record
    if (lines_since_phred_line): continue
    for phred_char in line.rstrip():
        nc_count+=1
        mean+=(math.pow(transform_base,ord(phred_char))-mean)/nc_count

print(mean/transform_scale)
