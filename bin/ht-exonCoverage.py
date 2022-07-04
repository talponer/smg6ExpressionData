#!/usr/bin/env python2.7

## This script is use for the analysis of mRNA stability. It counts
## reads in the whole gene locus (similar output as htseq-counts) or
## for reads that align to exon only. In both cases it output read
## counts at gene level.


import HTSeq
import collections
import argparse

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='This script calculate '
                                     'average read depth in exons.')
    parser.add_argument('-g','--gtf',
                        help='GTF file with genomic regions.',
                        nargs=1,
                        type=str)
    parser.add_argument('-b','--bam',
                        help='BAM file with aligned reads.',
                        nargs=1,
                        type=str)
    parser.add_argument('-c','--chr',
                        help='TAB delimited files with chromosome lengths',
                        nargs=1,
                        type=str)
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def matchReads(cigop, strd):
    global gene_ids, out, intervals, intervalStrand, readDepth
    if debug:
        print '  Cigar strand: ' + cigop.ref_iv.strand
        print '  Cigar type: ' + cigop.type
        print '  Cigar size:', cigop.size
        print '  Out is set to: ', out
    if strd == 'first':
        if cigop.ref_iv.strand == '+':
            cigop.ref_iv.strand = '-'
        else:
            cigop.ref_iv.strand = '+'
    if cigop.type != "M":
        if debug: print '  Skip this part'
        return None
    elif out:
        if debug: print '  This read was marked as out'
        pass
    else:
        for iv, val in intervals[ cigop.ref_iv ].steps():
            if len(val) == 0:
                out = 1
                if debug: print '  read matches outside region'
            for gid in val:
                if debug: print '  working on ' + gid
                if debug: print '    pair = ' + strd
                if cigop.ref_iv.strand == intervalStrand[ gid ]:
                    if debug: print '    the strand is good'
                    if not iv.contains(cigop.ref_iv):
                        if debug: print '    the reads span larger region'
                        counts[ "__alignment_partially_outside" ] += 1
                        out = 1
                        # if gid in gene_ids:
                        #     if debug: print '    the gene was already present, deleting it'
                        #     gene_ids.remove(gid)
                        continue
                    else:
                        if debug: print '    the read is contanined'
                        readDepth[cigop.ref_iv] += 1
                        # gene_ids.add(gid)
                else:
                    if debug: print '    strand is wrong: ' + first_almnt.iv.strand + \
                        ' ' + intervalStrand[ gid ]


def Average(lst):
    return sum(lst) / len(lst)

args = parserFunction()

debug   = args.debug
gtfFile = args.gtf[0]
bamFile = args.bam[0]
chrFile = args.chr[0]

fType = "exon"

## Read chr file:
chromlens = {}
chr = open(chrFile, 'r') 
Lines = chr.readlines() 
for line in Lines: 
    chrName, chrLength = line.strip().split('\t')
    chromlens[ chrName ] = int(chrLength)


## if debug: print fType, gtfFile, bamFile
    
## Read GTF file
gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
intervals = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
readDepth = HTSeq.GenomicArray( chromlens, stranded=False, typecode="i" )
counts = collections.Counter( )
intervalStrand = dict()

for feature in gtf_file:
    if feature.type == fType:
        intervals[ feature.iv ] += feature.attr["gene_id"]
        counts[ feature.attr["gene_id"] ] = 0
        intervalStrand[ feature.attr["gene_id"] ] = feature.iv.strand
        # intervals1[ feature.iv ] += feature.attr["gene_id"]
        # intervals2[ feature.iv ] += feature.attr["exon_id"]


## Read alignement file
almnt_file = HTSeq.BAM_Reader( bamFile )
for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=True ):
    if len(bundle) != 1:
        counts[ "__alignment_not_unique" ] += 1
        continue
    first_almnt, second_almnt = bundle[0]
    ## Check if read is mapped:
    if not first_almnt.aligned or not second_almnt.aligned:
        counts[ "__unmapped" ] += 1
        continue
    ## Check if quality is good:
    if first_almnt.aQual < 4 and second_almnt.aQual < 4:
        counts[ "__too_low_aQual" ] += 1
        continue
    
    ## Check where the read map accounting for possible gaps:
    gene_ids = set()
    out = 0    
    for cigop in first_almnt.cigar:
        pair = 'first'
        if debug: print 'Read ' + first_almnt.read.name + ' ' + pair + ' ' + first_almnt.iv.strand
        matchReads(cigop, pair)
    for cigop in second_almnt.cigar:
        pair = 'second'
        if debug: print 'Read ' + second_almnt.read.name + ' ' + pair + ' ' + second_almnt.iv.strand
        matchReads(cigop, pair)


## Count average read depth per exon:
for feature in gtf_file:
    if feature.type == fType:
        print feature.attr["exon_id"], feature.attr["gene_id"], Average(list(readDepth[ feature.iv ]))

    

