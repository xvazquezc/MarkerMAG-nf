#!/usr/bin/env python3

import os
import glob
import shutil
import argparse
from Bio import SeqIO
from datetime import datetime


polish_16s_usage = '''
===================== polish_16s example commands =====================

MarkerMAG polish_16s -in Matam_output.fa -out Matam_output_polished.fa

=======================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def polish_16s(args):

    file_in = args['in']
    file_out_ffn = args['out']

    file_out_path, file_out_base, file_out_ext = sep_path_basename_ext(file_out_ffn)

    barrnap_stdout   = '%s/%s.log'    % (file_out_path, file_out_base)
    file_out_gff     = '%s/%s.gff'    % (file_out_path, file_out_base)
    barrnap_cmd = 'barrnap --quiet %s 2> %s > %s' % (file_in, barrnap_stdout, file_out_gff)
    os.system(barrnap_cmd)

    input_seqs = SeqIO.to_dict(SeqIO.parse(file_in, 'fasta'))
    wrote_id = []
    file_out_ffn_handle = open(file_out_ffn, 'w')
    with open(file_out_gff) as gff_handle:
        for line in gff_handle:
            if line.startswith('#'):
                continue

            fields = line.rstrip().split('\t')
            if len(fields) != 9:
                continue

            seq_id, _, _, start, end, _, strand, _, attributes = fields
            if ('Name=16S_rRNA' not in attributes) and ('product=16S ribosomal RNA' not in attributes):
                continue

            each_16s = input_seqs[seq_id][int(start) - 1:int(end)]
            if strand == '-':
                each_16s = each_16s.reverse_complement()

            output_id = seq_id
            if seq_id in wrote_id:
                output_id = '%s_%s' % (seq_id, wrote_id.count(seq_id) + 1)

            file_out_ffn_handle.write('>%s\n' % output_id)
            file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
            wrote_id.append(seq_id)

    file_out_ffn_handle.close()

    if os.path.isfile('%s.fai' % file_in):
        os.remove('%s.fai' % file_in)


if __name__ == '__main__':

    matam_16s_parser = argparse.ArgumentParser()

    matam_16s_parser.add_argument('-in',  required=True, help='input 16S sequences')
    matam_16s_parser.add_argument('-out', required=True, help='output sequences')

    args = vars(matam_16s_parser.parse_args())
    polish_16s(args)
