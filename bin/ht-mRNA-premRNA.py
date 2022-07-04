#!/usr/bin/env python3

import HTSeq
import collections
import argparse
from math import ceil as ceil

def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='This script calculates '
                                     'average read depth for mRNA (exon'
                                     ' mapping reads) and pre-mRNA (intron'
                                     ' mapping reads). It output their '
                                     'averages at gene level in two files')
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
    parser.add_argument('-s','--select',
                        help='TAB delimited files list of selected exon and gene IDs.'
                        ' Default All.',
                        nargs=1,
                        type=str,
                        default='a')
    parser.add_argument('-m','--mrna',
                        help='File name for staring mRNA read counts.'
                        ' Default mRNAcounts',
                        nargs=1,
                        type=str,
                        default='mRNAcounts')
    parser.add_argument('-p','--premrna',
                        help='File name for staring pre-mRNA read counts. '
                        'Default premRNAcounts',
                        nargs=1,
                        type=str,
                        default='premRNAcounts')
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def matchReads(cigop, strd):
    global gene_ids, out, exonInter, intervalStrand, exonDepth, intron
    if debug:
        print('  Cigar strand:', cigop.ref_iv.strand)
        print('  Cigar type:', cigop.type)
        print('  Cigar size:', cigop.size)
        print('  Out is set to:', out)
    if strd == 'first':
        if cigop.ref_iv.strand == '+':
            cigop.ref_iv.strand = '-'
        else:
            cigop.ref_iv.strand = '+'
    if cigop.type != "M":
        if debug: print('  Skip this part')
        return None
    elif out:
        if debug: print('  This read was marked as out')
        pass
    else:
        for iv, gene in geneInter[ cigop.ref_iv ].steps():
            if len(gene) == 0: ## if read map outside gene stop
                out = 1
                if debug: print('  read matches outside gene')
                continue
            for gid in gene:
                if debug: print('  working on', gid)
                if debug: print('    pair =', strd)
                for eiv, egene in exonInter[ cigop.ref_iv ].steps():
                    if cigop.ref_iv.strand == intervalStrand[ gid ]:
                        mapRegion.append(cigop.ref_iv)
                        if debug:
                            print('    the strand is good')
                            print('    interval =', cigop.ref_iv)                    
                        if not eiv.contains(cigop.ref_iv): # if read does not map in exon
                            if debug: print('    read overlaps intron')
                            # intronDepth[ cigop.ref_iv ] += 1
                            intron = 1
                        else:
                            if debug: print('    read map exon')
                            # exonDepth[cigop.ref_iv] += 1
                            # gene_ids.add(gid)
                    else:
                        if debug: print('    strand is wrong:', first_almnt.iv.strand,
                                        intervalStrand[ gid ])


def Average(lst):
    if len(lst) > 0:
        return sum(lst) / len(lst)
    else:
        return 0

def printReadDepth(ftype, gtf, depth, fileName):
    f = open(fileName, 'w')
    if ftype == 'exon':
        nid = selExon
        tar = 'exon_id'
    else:
        nid = selGene
        tar = 'gene_id'
    oldGid = ''
    depthList = []
    totDepth = 0
    for feature in gtf:
        if feature.type == ftype:
            if (len(nid) > 0):
                if feature.attr[ tar ] not in nid:
                    continue
            gid = feature.attr["gene_id"]
            if gid != oldGid:
                value = str(ceil(Average(depthList)))
                string = oldGid + "\t" + value + "\n"
                if oldGid != '':
                    f.write(string)
                    if debug: print('END:', string, end='')
                oldGid = gid
                depthList = list(depth[ feature.iv ])
                totDepth = sum(list(depth[ feature.iv ]))
                if debug: print('BEG:', gid, ftype, str(ceil(Average(depthList))))
            else:
                depthList.append(Average(list(depth[ feature.iv ])))
                totDepth += sum(list(depth[ feature.iv ]))
                if debug: print('CON:', gid, ftype,
                                str(ceil(Average(list(depth[ feature.iv ])))),
                                str(ceil(Average(depthList))))
    string = gid + "\t" + str(ceil(Average(depthList))) + "\n"
    f.write(string)
    f.close()


def printCoverage(ftype, gInterval, depth, fileName):
    global exonInter
    f = open(fileName, 'w')
    ## count coverage for exon:
    if debug:
        print('Printing gene exon intervals:')
    for gid in sorted(gInterval):
        iv = gInterval[ gid ]
        depthList = []
        if debug: print(gid)
        if ftype == 'exon':
            for iv2, val in exonInter[ iv ].steps():
                if debug: print('  ', iv2, sorted(val))
                ## evaluate depth only if exon is for 1 gene only
                if len(val) == 1:
                    depthList.append(Average(list(depth[ iv2 ])))
        else:
            for iv2, val in exonInter[ iv ].steps():
                if debug: print('  ', iv2, sorted(val))
                ## evaluate depth only if exon is for 1 gene only
                if len(val) == 0:
                    depthList.append(Average(list(depth[ iv2 ])))
                    ## value = str(ceil(Average(depthList)))
            ## value = str(ceil(Average(list(depth[ iv ]))))
        value = str(ceil(Average(depthList)))
        string = gid + "\t" + value + "\n"
        f.write(string)
    f.close()


args = parserFunction()

debug   = args.debug
gtfFile = args.gtf[0]
bamFile = args.bam[0]
chrFile = args.chr[0]
selFile = args.select[0]
mrnaFile = args.mrna[0]
premrnaFile = args.premrna[0]

## Read chr file:
chromlens = {}
chr = open(chrFile, 'r') 
Lines = chr.readlines() 
for line in Lines: 
    chrName, chrLength = line.strip().split('\t')
    chromlens[ chrName ] = int(chrLength)

## Read exon list:
selExon = set()
selGene = set()
if (selFile != 'a'):
    if debug: print('Restrict counting')
    sel = open(selFile, 'r') 
    Lines = sel.readlines() 
    for line in Lines: 
        exonID, geneID = line.strip().split('\t')
        selExon.add(exonID)
        selGene.add(geneID)
    if debug:
        print('  got', len(selExon), 'unique exons')
        print('  got', len(selGene), 'unique genes')

## if debug: print(fType, gtfFile, bamFile
    
## Read GTF file
gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
exonInter = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
geneInter = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
exonDepth = HTSeq.GenomicArray( chromlens, stranded=False, typecode="i" )
intronDepth = HTSeq.GenomicArray( chromlens, stranded=False, typecode="i" )
counts = collections.Counter( )
intervalStrand = dict()
gInterval = dict()
eInterval = dict()

for feature in gtf_file:
    if feature.type == 'exon':
        if (len(selExon) > 0):
            if feature.attr['exon_id'] not in selExon:
                continue        
        exonInter[ feature.iv ] += feature.attr["gene_id"]
        ## eInterval[ feature.attr["exon_id"] ] = feature.iv
        exonDepth[ feature.iv ] = 0
        intervalStrand[ feature.attr["gene_id"] ] = feature.iv.strand
    elif feature.type == 'gene':
        if (len(selGene) > 0):
            if feature.attr['gene_id'] not in selGene:
                continue        
        geneInter[ feature.iv ] += feature.attr["gene_id"]
        gInterval[ feature.attr["gene_id"] ] = feature.iv
        intronDepth[ feature.iv ] = 0

        
## Read alignment file:
## almnt_file = HTSeq.BAM_Reader( bamFile )
almnt_file = HTSeq.SAM_Reader( bamFile )
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
    intron = 0
    mapRegion = list()
    for cigop in first_almnt.cigar:
        pair = 'first'
        if debug: print('Read', first_almnt.read.name, pair, first_almnt.iv.strand)
        matchReads(cigop, pair)
    for cigop in second_almnt.cigar:
        pair = 'second'
        if debug: print('Read', second_almnt.read.name, pair, second_almnt.iv.strand)
        matchReads(cigop, pair)
    if not out:
        if debug:
            print('Read', first_almnt.read.name, 'was marked as mapping a gene')
            ## print('Mapped regions:', mapRegion)
            if intron:
                print('  it is an intron read!')
            else:
                print('  it is an exon read!')
        for iv in mapRegion:
            if intron:
                intronDepth[iv] += 1
                if debug: print(iv, list(intronDepth[iv]))
            else:
                exonDepth[iv] += 1
                if debug: print(iv, list(exonDepth[iv]))
        if debug: print('')

## Count average read depth per exon:
## printReadDepth('exon', gtf_file, exonDepth, mrnaFile)
printCoverage('exon', gInterval, exonDepth, mrnaFile)

## printReadDepth('gene', gtf_file, intronDepth, premrnaFile)
printCoverage('gene', gInterval, intronDepth, premrnaFile)


    

