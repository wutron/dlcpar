
import os
import subprocess
import sys

data_dir = 'data/10e6-1x/'
log_file = '/research/yjw/jgardi/projects/dlcparIlp/examples/results.tsv'


def count_dups_and_coals(output):
    dups = 0
    coals = 0
    # I skip the first line because it has the column headers
    for line in output.split('\n')[1:]:
        columns = line.split('\t')
        dups += int(columns[4])
        coals += int(columns[6])

    return dups, coals


for sub_dir_name in sorted([int(dir) for dir in os.listdir(data_dir) if not dir.startswith('.')]):
    print(sub_dir_name)
    coal_tree = '{data_dir}{sub_dir_name}/{sub_dir_name}.coal.tree'.format(**vars())

    dlcpar_command_base = ['dlcpar', '-s', 'config/flies.stree', '-S', 'config/flies.smap', '-I', '.coal.tree',
                           '-O', '.dlcpar', '-x1234', coal_tree]
    try:
        subprocess.run(dlcpar_command_base + ['-K 0.001'])
    except Exception:
        continue

    results_command = 'echo ' +  coal_tree + ' | tree-events-dlc -s config/flies.stree '\
                                            '-S config/flies.smap -T .coal.tree'
    output_with_low_K = subprocess.getoutput(results_command)
    subprocess.run(dlcpar_command_base +  ['-K 10'])
    output_with_high_K = subprocess.getoutput(results_command)

    dups_with_low_K, coals_with_low_K = count_dups_and_coals(output_with_low_K)
    dups_with_high_K, coals_with_high_K = count_dups_and_coals(output_with_high_K)

    if dups_with_low_K != dups_with_high_K:
        print('dups are different', dups_with_low_K, dups_with_high_K)

    if coals_with_low_K != coals_with_high_K:
        print('coals are different', coals_with_low_K, coals_with_high_K)



