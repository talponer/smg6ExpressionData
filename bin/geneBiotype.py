#!/usr/bin/env python3

import HTSeq
import collections
import argparse
from math import ceil as ceil

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='This script check if a gene '
                                     'belong to the following chategories:'
                                     '1. NMD (if one isoform is NMD sensitive);'
                                     '2. Retain intron (if one isoform is RI);'
                                     '3. Protein coding if at least one isoform'
                                     ' is so')
    parser.add_argument('-g','--gtf',
                        help='GTF file with genomic regions.',
                        nargs=1,
                        type=str)
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def parseGtf(gtfFile, debug):
    geneBioType = dict()
    nmd = 0
    ri = 0
    gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
    for feature in gtf_file:
        if feature.type == 'transcript':
            gID = feature.attr['gene_id']
            tBioType = feature.attr['transcript_biotype']
            if debug: print('Gene ID:', gID)
            if gID not in geneBioType:
                geneBioType[gID] = tBioType
                if debug: print('  biotype:', tBioType)
            else:
                if tBioType == "nonsense_mediated_decay":
                    geneBioType[gID] = tBioType
                    nmd += 1
                    if debug: print('  biotype:', tBioType)
                elif tBioType == "retained_intron":
                    if geneBioType[gID] != "nonsense_mediated_decay":
                        geneBioType[gID] = tBioType
                        ri += 1
                        if debug: print('  biotype:', tBioType)
                elif tBioType == "protein_coding":
                    if debug: print('Existing biotype:', geneBioType[gID])
                    if geneBioType[gID] == "nonsense_mediated_decay" or geneBioType[gID] == "retained_intron":
                        continue
                    else:
                        geneBioType[gID] = tBioType
                        if debug: print('  changing to:', tBioType)
    return(geneBioType, nmd, ri)

args = parserFunction()

debug   = args.debug
gtfFile = args.gtf[0]
geneBioType, nmd, ri = parseGtf(gtfFile, debug)


## print('NMD:', nmd)
## print('RI:', ri)
print("GeneID", "biotype", sep='\t') 
for id in sorted(geneBioType):
    print(id, geneBioType[id], sep='\t')
