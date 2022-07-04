#!/usr/bin/env python3

## This script is use for the analysis of mRNA stability. It counts
## reads in the whole gene locus (similar output as htseq-counts) or
## for reads that align to exon only. In both cases it output read
## counts at gene level.


import HTSeq
import collections
import argparse

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='This script counts reads '
                                     'that map to gene locuses (exon+intron) '
                                     'or only to exons. The output is '
                                     'always summariezed at the gene level.'
                                     'Alignemnts must be in BAM format.')
    parser.add_argument('bed',
                        help='BED file with genomic regions.',
                        nargs=1,
                        type=str)
    parser.add_argument('bam',
                        help='BAM file with aligned reads.',
                        nargs=1,
                        type=str)
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def matchReads(cigop, strd):
    global gene_ids, out, intervals, intervalStrand
    if debug:
        print('  Cigar strand: ' + cigop.ref_iv.strand)
        print('  Cigar type: ' + cigop.type)
        print('  Cigar size:', cigop.size)
        print('  Out is set to: ', out)
    if strd == 'first':
        if cigop.ref_iv.strand == '+':
            cigop.ref_iv.strand = '-'
        else:
            cigop.ref_iv.strand = '+'
    if cigop.type != "M":
        if debug: print('  Skip this part')
        return None
    else:
        good = 0
        gid = ''
        for iv, val in intervals[ cigop.ref_iv ].steps():
            for gid in val:
                if debug: print('  working on ' + gid)
                if debug: print('    pair = ' + strd)
                if cigop.ref_iv.strand == intervalStrand[ gid ]:
                    if debug: print('    the strand is good')
                    if debug: print('    feature interval:', iv)
                    if debug: print('    read interval:', cigop.ref_iv)
                    if debug: print('    the read is contanined')
                    gene_ids.add(gid)
                    break
                else:
                    if debug: print('    strand is wrong: ' + first_almnt.iv.strand + \
                        ' ' + intervalStrand[ gid ])


isExon = 0
args = parserFunction()

debug   = args.debug
bedFile = args.bed[0]
bamFile = args.bam[0]

## Read GTF file
intervals = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
counts = collections.Counter( )
intervalStrand = dict()
totGene = set()

if debug: print('Reading BED file')

for line in open(bedFile):
    l = line.strip('\n')
    if debug: print('  L:', l)
    fields = l.split( "\t" )
    iv = HTSeq.GenomicInterval( fields[0], int(fields[1]), int(fields[2]) )
    intervals[ iv ] += fields[3]
    intervalStrand[ fields[3] ] = fields[5]
    counts[ fields[3] ] = 0
    

if debug: print('Reading BAm file')
almnt_file = HTSeq.BAM_Reader( bamFile )
for bundle in HTSeq.pair_SAM_alignments( almnt_file, bundle=False ):
    if debug: print('  bundle:', bundle)
    
    first_almnt, second_almnt = bundle
    if debug:
        print('  first:', first_almnt)
        print('  second:', second_almnt)

        
    ## Check where the read map accounting for possible gaps:
    gene_ids = set()
    out = 0
    if first_almnt:
        if first_almnt.not_primary_alignment:
            if debug: print('  FA.PA:', first_almnt.not_primary_alignment)
            continue
        for cigop in first_almnt.cigar:
            pair = 'first'
            if debug: print('Read ' + first_almnt.read.name + ' ' + pair + ' ' + first_almnt.iv.strand)
            matchReads(cigop, pair)
    if second_almnt:
        if second_almnt.not_primary_alignment:
            if debug: print('  SA.PA:', second_almnt.not_primary_alignment)
            continue
        for cigop in second_almnt.cigar:
            pair = 'second'
            if debug: print('Read ' + second_almnt.read.name + ' ' + pair + ' ' + second_almnt.iv.strand)
            matchReads(cigop, pair)
    ## Make sure to not count reads that map outside regions in any of
    ## the pairs
    if out: gene_ids = set()
    
    ## Count the reads
    if len(gene_ids) == 1:
        gene_id = list(gene_ids)[0]
        counts[ gene_id ] += 1
    # elif len(gene_ids) == 0:
    #     counts[ "__no_feature" ] += 1
    # else:
    #     counts[ "__ambiguous" ] += 1

for gene_id in sorted(counts.items()):
    print(gene_id[0], gene_id[1], sep="\t")
