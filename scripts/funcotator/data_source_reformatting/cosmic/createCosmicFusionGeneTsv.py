#!/usr/bin/env python

from shared_utils import count_lines

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import csv
from oncotator.utils.GenericTsvReader import GenericTsvReader 
from collections import OrderedDict
import re
def parseOptions():
    
    # Process arguments
    desc = ''' Create the COSMIC fusion gene Oncotator datasource.'''
    epilog = ''' NOTE: This script will load large portions of the COSMIC input file into RAM.

    Updated for v76
    '''
    parser = ArgumentParser(description=desc, formatter_class=RawDescriptionHelpFormatter, epilog=epilog)
    parser.add_argument("ds_file", type=str, help="COSMIC datasource filename. For example, 'CosmicFusionExport.tsv' ")
    parser.add_argument("output_file", type=str, help="TSV filename for output.  File will be overwritten if it already exists.")
    
    args = parser.parse_args()
    return args

def renderFusionGeneDictEntry(geneKey, fusionGeneDict):
    fusionGeneSubDict = fusionGeneDict[geneKey]
    resultList = []
    for k in fusionGeneSubDict.keys():
        print k
        print str(fusionGeneSubDict[k])
        summaryString = "%(k)s(%(fgene)s)" % {'k':k,'fgene':str(fusionGeneSubDict[k])}
        resultList.append(summaryString)
    return '|'.join(resultList)

if __name__ == '__main__':
    args = parseOptions()
    inputFilename = args.ds_file
    outputFilename = args.output_file
    num_lines = count_lines(inputFilename)
    tsvReader = GenericTsvReader(inputFilename)
    headers = tsvReader.getFieldNames()
    print('Found headers (input): ' + str(headers))
    if "Translocation Name" not in headers:
        raise NotImplementedError("Could not find Translocation Name column in the input file.")
    
    outputHeaders = ['gene', 'fusion_genes', 'fusion_id']
    
    # Create a dictionary where key is the gene and value is another dict: {fusion_gene:count}
    fusionGeneDict = OrderedDict()
    last_i = 0
    for i, line in enumerate(tsvReader):
        fusion_gene_description = line['Translocation Name']
        
        if len(fusion_gene_description.strip()) == 0:
            # blank
            continue

        # geneListKeys = fusionGene.split('/')

        genes_in_this_fusion = re.findall("_*([A-Z0-9\-\.]+)\{", fusion_gene_description)

        for k in genes_in_this_fusion:
            if k not in fusionGeneDict.keys():
                fusionGeneDict[k] = dict()

            # Look for the fusion
            if fusion_gene_description not in fusionGeneDict[k].keys():
                fusionGeneDict[k][fusion_gene_description] = 0

            fusionGeneDict[k][fusion_gene_description] = fusionGeneDict[k][fusion_gene_description] + 1

        if i - last_i > round(float(num_lines)/100.0):
            print("{:.0f}% complete".format(100 * float(i)/float(num_lines)))
            last_i = i
        
    # Render the fusionGeneDict
    tsvWriter = csv.DictWriter(file(outputFilename,'w'), outputHeaders, delimiter='\t', lineterminator="\n")
    tsvWriter.fieldnames = outputHeaders
    tsvWriter.writeheader()
    for k in fusionGeneDict.keys():
        
        row = dict()
        row['gene'] = k
        row['fusion_genes'] = renderFusionGeneDictEntry(k, fusionGeneDict) 
        tsvWriter.writerow(row)
    
    pass
    
    
