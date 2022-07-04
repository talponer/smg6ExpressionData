#!/usr/bin/env python3

import HTSeq
import collections
import argparse
from math import ceil as ceil

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='This script calculates '
                                     'the mRNA and pre-mRNA (intron) lengths'
                                     ' for each gene.')
    parser.add_argument('-g','--gtf',
                        help='GTF file with genomic regions.',
                        nargs=1,
                        type=str)
    parser.add_argument('-l','--lfile',
                        help='File with expressed exons. This is a two '
                        'columns file with exon ID and genes ID',
                        nargs=1,
                        type=str)
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def parseGtf(gtfFile, exons, debug):
    mrnaLength = dict()
    premrnaLength = dict()
    geneLength = dict()
    intervals = dict()
    geneStart = dict()
    geneEnd = dict()
    geneChr = dict()
    gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
    for feature in gtf_file:
        if feature.type == 'exon':
            gID = feature.attr['gene_id']
            eID = feature.attr['exon_id']
            if gID in exons:
                if eID not in exons[gID]:
                    continue
                else:
                    interval = feature.iv
                    if gID not in intervals:
                        chr = list(interval.chrom)
                        intervals[gID] = HTSeq.GenomicArrayOfSets( 'auto',
                                                                  stranded=False )
                        intervals[gID][interval] += gID
                    else:
                        intervals[gID][interval] += gID
        elif feature.type == 'gene':
            gID = feature.attr['gene_id']
            if gID in exons:
                eInt = feature.iv.length
                geneLength[gID] = feature.iv.length
                geneStart[gID] = feature.iv.start
                geneEnd[gID] = feature.iv.end
                geneChr[gID] = feature.iv.chrom
                # if debug:
                #     print('Gene ID:', gID)
                #     print('Length :', geneLength[gID])
    for ID in sorted(intervals):
        if debug: print(ID)
        mrnaLength[ID] = 0
        gas = intervals[ID]
        region = HTSeq.GenomicInterval( geneChr[ID], geneStart[ID],
                                        geneEnd[ID], "." ) 
        for iv, value in gas[ region ].steps():
            if len(value) > 0:
                if list(value)[0] == ID:
                    mrnaLength[ID] += iv.length
                    if debug:
                        print('Added :', iv.length)
    return(mrnaLength, geneLength)


def overlap(iv1, iv2):
    s1 = iv1.start
    e1 = iv1.end
    s2 = iv2.start
    e2 = iv2.end
    if(s1 > s2):
        eInt = e1 - e2
    else:
        eInt = s2 - s1
    return(eInt)

def getExons(lfile, debug):
    exons = dict()
    lst = []
    if debug: print('Restrict counting')
    sel = open(lfile, 'r') 
    Lines = sel.readlines() 
    for line in Lines: 
        exonID, geneID = line.strip().split('\t')
        if geneID in exons:
            exons[geneID].add(exonID)
        else:
            exons[geneID] = set(exonID)
    return exons
        

args = parserFunction()

debug   = args.debug
gtfFile = args.gtf[0]
listFile = args.lfile[0]

exons = getExons(listFile, debug)
mLength, gLength = parseGtf(gtfFile, exons, debug)

print("GeneID", "ExonLength", "IntronLength", sep='\t') 
for id in sorted(mLength):
    print(id, mLength[id], gLength[id] - mLength[id], sep='\t')
