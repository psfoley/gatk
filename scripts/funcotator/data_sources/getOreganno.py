#!/usr/bin/env bash

########################################################################
# Imports:

import csv
import urllib
from GenericTsvReader import GenericTsvReader

########################################################################
# Constants:

FILE_URL = 'http://www.oreganno.org/dump/ORegAnno_Combined_2016.01.19.tsv'
RAW_FILE_NAME = FILE_URL.split('/')[-1]

OUTPUT_HEADERS = ['Build', 'Chr', 'Start', 'End', 'ID', 'Values']

OUT_HG19_FILE_NAME_TEMPLATE = 'oreganno_%s.HG19.tsv'
OUT_HG38_FILE_NAME_TEMPLATE = 'oreganno_%s.HG38.tsv'

VALUES_DELIMITER = '|'

########################################################################
# Functions:


def getValuesDataFromRowDict(r):
    values = ''

    for key in ['Outcome', 'Type', 'Gene_Symbol', 'Gene_ID', 'Gene_Source', 'Regulatory_Element_Symbol',
                'Regulatory_Element_ID', 'Regulatory_Element_Source', 'dbSNP_ID', 'PMID', 'Dataset']:
        values = values + key + "=" + r[key]

    return values

########################################################################
# Main:

if __name__ == '__main__':

    # Download the Raw file:
    urllib.urlretrieve(FILE_URL, RAW_FILE_NAME)

    # Get the version of the file:
    file_version = RAW_FILE_NAME.replace('ORegAnno_Combined_', '').replace('.tsv', '').replace('.', '')

    # Set up the output files:
    out_hg19_tsv_writer = csv.DictWriter(file(OUT_HG19_FILE_NAME_TEMPLATE % file_version, 'w'), OUTPUT_HEADERS, delimiter='\t', lineterminator="\n")
    out_hg19_tsv_writer.fieldnames = OUTPUT_HEADERS
    out_hg19_tsv_writer.writeheader()
    out_hg38_tsv_writer = csv.DictWriter(file(OUT_HG38_FILE_NAME_TEMPLATE % file_version, 'w'), OUTPUT_HEADERS, delimiter='\t', lineterminator="\n")
    out_hg38_tsv_writer.fieldnames = OUTPUT_HEADERS
    out_hg38_tsv_writer.writeheader()

    # Now that we have the file, go through it and reformat it:
    tsvReader = GenericTsvReader(RAW_FILE_NAME)
    headers = tsvReader.getFieldNames()
    print('Found headers (input): ' + str(headers))

    # Go through our file:
    for i, line in enumerate(tsvReader):

        if line['Species'].lower() != 'Homo sapiens':
            print 'Ignoring non-human record on line:', i
            continue

        # Get the trivial fields here:
        row = dict()
        row['Build'] = line['Build']
        row['Chr'] = line['Chr']
        row['Start'] = line['Start']
        row['End'] = line['End']
        row['ID'] = line['ORegAnno_ID']

        # Get the values field from our helper method:
        row['Values'] = getValuesDataFromRowDict(row)

        if line['Build'].lower() == 'hg19':
            out_hg19_tsv_writer.writerow(row)
        elif line['Build'].lower() == 'hg38':
            out_hg38_tsv_writer.writerow(row)
        else:
            print 'Ignoring incompatible build record on line:', i



