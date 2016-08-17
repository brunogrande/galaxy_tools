#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import shlex
import time


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--tumor_id')
    parser.add_argument('--normal_id')
    parser.add_argument('--input_vep')
    parser.add_argument('--output_maf')

    args = parser.parse_args()
    args_dict = vars(args)
    cwd = os.getcwd()
    exit_codes = []
    print os.environ

    # Since vcf2maf.pl creates an intermediate file in the same directory as the source
    # VCF file, we need to link to the VCF file in the current working directory.
    args_dict['link_to_input_vep'] = cwd + '/input.vep.vcf'
    command_0 = ('ln -s {input_vep} {link_to_input_vep}').format(**args_dict)
    print command_0
    command_0_ready = shlex.split(command_0)
    process_0 = subprocess.Popen(command_0_ready)
    process_0_status = process_0.wait()
    exit_codes.append(process_0_status)

    # Running vcf2maf.pl
    command_1 = ('/home/bgrande/galaxy-dist/tools/gtools/vcf2maf/vcf2maf.pl '
                 '--input-vep {link_to_input_vep} '
                 '--output-maf {output_maf} '
                 '--vep-path /home/bgrande/software/variant_effect_predictor '
                 '--vep-data /home/bgrande/.vep '
                 '--tumor-id {tumor_id} '
                 '--normal-id {normal_id}'
                 ).format(**args_dict)
    print command_1
    command_1_ready = shlex.split(command_1)
    process_1 = subprocess.Popen(command_1_ready)
    process_1_status = process_1.wait()
    exit_codes.append(process_1_status)

    # Only exit with status code 0 if all subprocesses exited with a 0 status
    if any(exit_codes):  # If any exit codes are not 0 (i.e., False)
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == '__main__':
    main()
