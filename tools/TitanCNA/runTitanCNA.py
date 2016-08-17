#!/usr/bin/env python


import argparse
import subprocess
import os
import re
import sys


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--tumour_bam')
    parser.add_argument('--tumour_name')
    parser.add_argument('--normal_bam')
    parser.add_argument('--normal_name')
    parser.add_argument('--output')
    args = parser.parse_args()
    args_dict = vars(args)
    args_dict['normal_name'] = 'normal'
    args_dict['tumour_name'] = 'tumour'
    number_of_clusters = 2
    exit_codes = []

    # Find TITANRunner's directory
    command_0 = 'which TitanRunner.py'
    command_0_output = subprocess.check_output(command_0, shell=True)
    args_dict['titanrunner_dir'] = os.path.dirname(command_0_output)

    # While debugging, use previously generated files
    # os.chdir('/home/bgrande/results/galaxy_output')

    # Command 1
    # Identify heterozygous germline positions (and their genotypes) in normal BAM file
    command_1 = ('samtools view -b -F 0x0400 {normal_bam} | '
                 'samtools mpileup -q 20 -Q 10 -uf /home/bgrande/data/homo_sapiens_hg19.chr22.fa - | '
                 'bcftools view -cgv - | '
                 'python {titanrunner_dir}/scripts/apollohFindHetPsns.py normal_hetpsns.txt;'
                 ).format(**args_dict)
    print command_1
    commmand_1_status = subprocess.call(command_1, shell=True)
    exit_codes.append(commmand_1_status)

    # Command 2
    # Calculate read counts for both normal and tumour BAM files
    command_2 = ('{titanrunner_dir}/scripts/readCounter -w 1000 -q 0 {normal_bam} | '
                 'python {titanrunner_dir}/scripts/filter_chromosomes.py > '
                 'normal_read_counts.wig; '
                 '{titanrunner_dir}/scripts/readCounter -w 1000 -q 0 {tumour_bam} | '
                 'python {titanrunner_dir}/scripts/filter_chromosomes.py > '
                 'tumour_read_counts.wig'
                 ).format(**args_dict)
    print command_2
    command_2_status = subprocess.call(command_2, shell=True)
    exit_codes.append(command_2_status)

    # Command 3
    # Correct read counts for GC content and mappability biases, and generate log ratios
    command_3 = ('R --no-save --args project_name '
                 'tumour_read_counts.wig normal_read_counts.wig '
                 '/home/bgrande/software/src/HMMcopy/data/gc_hg19.chr22.wig '
                 '/home/bgrande/software/src/HMMcopy/data/map_hg19.chr22.wig cn_data.txt < '
                 '{titanrunner_dir}/scripts/correctReads.R 2>&1'
                 ).format(**args_dict)
    print command_3
    command_3_status = subprocess.call(command_3, shell=True)
    exit_codes.append(command_3_status)

    # Command 4
    # Determine read counts for each germline heterozygous position identified previously
    command_4 = ('python {titanrunner_dir}/scripts/count.py normal_hetpsns.txt {tumour_bam} '
                 '/home/bgrande/data/homo_sapiens_hg19.chr22.fa 10 20 > '
                 'tumour_read_counts_at_hetpsns.txt'
                 ).format(**args_dict)
    print command_4
    command_4_status = subprocess.call(command_4, shell=True)
    exit_codes.append(command_4_status)

    # Command 5
    # Filter out positions that are not in dbSNP
    command_5 = ('python {titanrunner_dir}/scripts/filterCounts.py '
                 '-t tumour_read_counts_at_hetpsns.txt '
                 '-r /home/bgrande/data/common_all.vcf.gz '
                 '-o tumour_read_counts_at_hetpsns.filtered.txt'
                 ).format(**args_dict)
    print command_5
    command_5_status = subprocess.call(command_5, shell=True)
    exit_codes.append(command_5_status)

    # Command 6
    # Run TitanCNA successively with different number of clusters (1-5)
    for i in range(1, number_of_clusters + 1):
        command_6 = ('R --no-save --args {tumour_name}_vs_{normal_name} '
                     'tumour_read_counts_at_hetpsns.filtered.txt cn_data.txt '
                     '/home/bgrande/software/src/HMMcopy/data/map_hg19.chr22.wig '
                     '{0} 4 2 '
                     'titan_cluster_{0}.txt titan_cluster_{0}_params.txt '
                     '0 TRUE 0.5 map 50 1e-300 1e9 1e9 15000 15000 8 TRUE <'
                     '{titanrunner_dir}/scripts/titan.R'
                     ).format(i, **args_dict)
        print command_6
        command_6_status = subprocess.call(command_6, shell=True)
        exit_codes.append(command_6_status)

    # Command 7
    # Identify the optimal number of clusters (i.e., lowest S_Dbw index)
    indices = {}
    params_file_template = 'titan_cluster_{}_params.txt'
    for i in range(1, number_of_clusters + 1):
        params_file = open(params_file_template.format(i))
        for line in params_file:
            if line.startswith('S_Dbw validity index'):
                index_value = float(re.search(r'\d+\.\d+', line).group())
                indices[i] = index_value
    optimal_cluster_number = sorted(indices.items(), key=lambda x: x[1])[0][0]

    # Command 8
    # Create SEG files for the optimal number of clusters
    command_8 = ('{titanrunner_dir}/scripts/createTITANsegmentfiles.pl '
                 '-id={tumour_name}_vs_{normal_name} '
                 '-infile=titan_cluster_{0}.txt '
                 '-outfile={tumour_name}_vs_{normal_name}_segs.txt '
                 '-outIGV={tumour_name}_vs_{normal_name}.seg;'
                 ).format(optimal_cluster_number, **args_dict)
    print command_8
    command_8_status = subprocess.call(command_8, shell=True)
    exit_codes.append(command_8_status)

    # Command 9
    # Plot the CNAs
    command_9 = ('R --no-save --args '
                 '{tumour_name}_vs_{normal_name}_cluster_{0}.RData '
                 '. {0} {tumour_name}_vs_{normal_name} {tumour_name}_vs_{normal_name} < '
                 '{titanrunner_dir}/scripts/plot_titan.R'
                 ).format(optimal_cluster_number, **args_dict)
    print command_9
    command_9_status = subprocess.call(command_9, shell=True)
    exit_codes.append(command_9_status)

    # Command 10
    # Return plot PNG file as tool output
    command_10 = ('mv {tumour_name}_vs_{normal_name}_'
                  '{tumour_name}_vs_{normal_name}_cluster_{0}/'
                  '{tumour_name}_vs_{normal_name}_'
                  '{tumour_name}_vs_{normal_name}_cluster_{0}_chr22.png '
                  '{output}'
                  ).format(optimal_cluster_number, **args_dict)
    print command_10
    command_10_status = subprocess.call(command_10, shell=True)
    exit_codes.append(command_10_status)

    # Obtain snapshot of working directory
    # subprocess.call('rm /home/bgrande/results/galaxy_output/*; cp -R . /home/bgrande/results/galaxy_output/',
    #                 shell=True)

    # Only exit with status code 0 if all subprocesses exited with a 0 status
    if any(exit_codes):  # If any exit codes are not 0 (i.e., False)
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == '__main__':
    main()
