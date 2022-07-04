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
    parser = argparse.ArgumentParser(description='This script counts reads '
                                     'that map to gene locuses (exon+intron) '
                                     'or only to exons. The output is '
                                     'always summariezed at the gene level.'
                                     'Alignemnts must be in BAM format.')
    parser.add_argument('-e','--exon',
                        help='Count reads that map to exon only. Otherwise '
                        'count reads that map in the whole gene locus '
                        '(i.e. exon AND introns). Output at gene level',
                        action="store_true")
    parser.add_argument('-g','--gtf',
                        help='GTF file with genomic regions.',
                        nargs=1,
                        type=str)
    parser.add_argument('-b','--bam',
                        help='BAM file with aligned reads.',
                        nargs=1,
                        type=str)
    parser.add_argument('-s','--select',
                        help='TAB delimited files list of selected exon and gene IDs.'
                        ' Default All.',
                        nargs=1,
                        type=str,
                        default='a')
    parser.add_argument('-d','--debug',
                        help='Print verbose output for debugging',
                        action="store_true")
    return parser.parse_args()


def matchReads(cigop, strd):
    global gene_ids, out, intervals, intervalStrand
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
        good = 0
        gid = ''
        for iv, val in intervals[ cigop.ref_iv ].steps():
            if out and not good:
                if debug: print '  read was marked as out'
                break
            if len(val) == 0:
                out = 1
                good = 0
                if debug: print '  read matches outside region'
                if gid in gene_ids:
                    if debug: print '    read mapped also to a gene, deleting it'
                    gene_ids.remove(gid)
                break
            for gid in val:
                if debug: print '  working on ' + gid
                if debug: print '    pair = ' + strd
                if cigop.ref_iv.strand == intervalStrand[ gid ]:
                    if debug: print '    the strand is good'
                    if debug: print '    feature interval:', iv
                    if debug: print '    read interval:', cigop.ref_iv
                    if not iv.contains(cigop.ref_iv):
                        if debug: print '    the reads span larger region'
                        counts[ "__alignment_partially_outside" ] += 1
                        if not good:
                            out = 1
                            if gid in gene_ids:
                                if debug: print '    the gene was already present, deleting it'
                                gene_ids.remove(gid)
                        # continue
                    else:
                        if debug: print '    the read is contanined'
                        gene_ids.add(gid)
                        out = 0
                        good = 1
                        break
                else:
                    if debug: print '    strand is wrong: ' + first_almnt.iv.strand + \
                        ' ' + intervalStrand[ gid ]


isExon = 0
args = parserFunction()

isExon  = args.exon
debug   = args.debug
gtfFile = args.gtf[0]
bamFile = args.bam[0]
selFile = args.select[0]

## Read chr file:
selExon = set()
selGene = set()
if (selFile != 'a'):
    if debug: print 'Restrict counting'
    sel = open(selFile, 'r') 
    Lines = sel.readlines() 
    for line in Lines: 
        exonID, geneID = line.strip().split('\t')
        selExon.add(exonID)
        selGene.add(geneID)
    if debug:
        print '  got', len(selExon), 'unique exons'
        print '  got', len(selGene), 'unique genes'


if isExon:
    fType = "exon"
    selected = selExon
    feat = 'exon_id'
else:
    fType = "gene"
    selected = selGene
    feat = 'gene_id'

## if debug: print fType, gtfFile, bamFile
    
## Read GTF file
gtf_file = HTSeq.GFF_Reader( gtfFile, end_included=True )
intervals = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
counts = collections.Counter( )
intervalStrand = dict()
totGene = set()

if debug: print 'Reading GTF file'

for feature in gtf_file:
    if feature.type == fType:
        if (len(selected) > 0):
            if feature.attr[feat] not in selected:
                continue
        intervals[ feature.iv ] += feature.attr["gene_id"]
        counts[ feature.attr["gene_id"] ] = 0
        intervalStrand[ feature.attr["gene_id"] ] = feature.iv.strand
        totGene.add( feature.attr["gene_id"] )
        # intervals1[ feature.iv ] += feature.attr["gene_id"]
        # intervals2[ feature.iv ] += feature.attr["exon_id"]

if debug: print '  counting on', len(totGene), 'genes'

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
    ## Make sure to not count reads that map outside regions in any of
    ## the pairs
    if out: gene_ids = set()
    
    ## Count the reads
    if len(gene_ids) == 1:
        gene_id = list(gene_ids)[0]
        counts[ gene_id ] += 1
    elif len(gene_ids) == 0:
        counts[ "__no_feature" ] += 1
    else:
        counts[ "__ambiguous" ] += 1

for gene_id in sorted(counts.items()):
    print gene_id[0], gene_id[1]
