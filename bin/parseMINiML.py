#!/usr/bin/env python3

import xml.etree.ElementTree as ET

tree = ET.parse('/mnt/data/rdreos/Projects/ribuorf/publicData/takahashi12/GSE39977_family.xml')
root = tree.getroot()

for sample in root.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample'):
    sampleChip = ''
    sampleTime = ''
    sampleId = sample.get('iid')
    #print(sampleId)
    channel = sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Channel')
    sampleTitle = sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Title').text
    for char in channel.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics'):
        if (char.get('tag') == 'time'):
            sampleTime = char.text.replace("\n",'');
        if (char.get('tag') == 'chip antibody'):
            sampleChip = char.text.replace("\n",'');
    for sup in sample.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Supplementary-Data'):
        if(sup.get('type') == 'BED'):
            sampleFile = sup.text.replace("\n",'');
    print(sampleId, sampleTitle, sampleChip, sampleTime, sampleFile, sep='\t')
        
