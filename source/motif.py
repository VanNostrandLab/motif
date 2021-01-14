#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Identify motifs for CLIP peak clusters.
"""

import os
import subprocess
import sys
import argparse

ANNOTATE_PEAK = 'annotate_peaks_bedformat_with_proxdistal_lncRNA.pl'
RUN_HOMER = 'pull_seqs_for_regions_and_run_homer_clip_analysis_version.pl'


def main():
    parser = argparse.ArgumentParser(description=__doc__, prog='motif',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('BED', type=str,
                        help='Path to the BED file, usually this is the peak clusters normalized compressed bed file.')
    parser.add_argument('SPECIES', type=str,
                        help="Species shortname, can be one of these: hg19, hg39, mm10.")
    parser.add_argument('OUTDIR', type=str,
                        help="Path to the output directory, default: the BED file's parent directory.")
    parser.add_argument('UID', type=str,
                        help='User specified identifier for the dataset.')
    
    args = parser.parse_args()
    outdir = args.OUTDIR
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    bed = os.path.abspath(args.BED)
    annotated_bed = os.path.join(outdir, f'{args.UID}.peak.clusters.annotated.with.proxdist.miRlncRNA')
    if os.path.isfile(annotated_bed) and os.path.getsize(annotated_bed):
        print('Annotated peak cluster file already exists, skip re-annotate.')
    else:
        print('Annotating peak clusters ...')
        cmd = ['perl', ANNOTATE_PEAK, bed, args.SPECIES, annotated_bed]
        p = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stdout, universal_newlines=True)
        if p.returncode:
            sys.exit(p.returncode)

    print('Pulling sequences for regions and running homer ...')
    cmd = ['perl', RUN_HOMER, annotated_bed, args.SPECIES, args.UID]
    p = subprocess.run(cmd, stdout=sys.stdout, stderr=sys.stdout, universal_newlines=True)
    if p.returncode:
        sys.exit(p.returncode)
    

if __name__ == '__main__':
    main()
