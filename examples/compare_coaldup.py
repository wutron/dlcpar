
import os
import subprocess
import sys
import re
from glob import *

sim_flies_dir = '/research/yjw/jgardi/work/dlcpar_ilp/sim-flies'
data_sets = ['10e6-1x', '10e6-2x', '10e6-4x']
examples_dir = '/research/yjw/jgardi/projects/dlcparIlp/examples'

#==========
# helper functions

def get_runtime_and_cost(file):
    lines = file.readlines()
    runtime_line = [x for x in lines if 'Runtime:' in x][0]
    cost_line = [x for x in lines if 'Cost:' in x][0]

    return re.search(r'\d+\.\d+', runtime_line).group(), re.search(r'\d+\.\d+', cost_line).group()


def sorted_fams_in_data_set(data_set):
    """TODO"""
    sorted_fam_ids = sorted([int(d) for d in os.listdir(sim_flies_dir + '/' + data_set) if d.isdigit()])
    return [(data_set, str(fam_id)) for fam_id in sorted_fam_ids]


def test(fam_files_common_prefix):
    """
    :param str fam_files_common_prefix: the prefix sahred by all the files for this family
    :return tuple: the costs and runtime from testing on this family
    """
    coal_tree = fam_files_common_prefix + '.coal.tree'
    args = ['-s', examples_dir + '/config/flies.stree', '-S', examples_dir + '/config/flies.smap', coal_tree]
    # dlcpar_command_base = ['dlcpar'] + args + ['-I',  '.coal.tree']

    # subprocess.run(dlcpar_command_base + ['-K 0.001'])
    # dlcpar_info_file = open(fam_files_common_prefix + '.dlcpar.info', 'r')
    # low_k_runtime, low_k_cost = get_runtime_and_cost(dlcpar_info_file)
    # dlcpar_info_file.close()

    # subprocess.run(dlcpar_command_base + ['-K 10'])
    # dlcpar_info_file = open(fam_files_common_prefix + '.dlcpar.info', 'r')
    # high_k_runtime, high_k_cost = get_runtime_and_cost(dlcpar_info_file)
    # dlcpar_info_file.close()

    # if high_k_cost != low_k_cost:
    #     print(('for {fam_files_common_prefix} the cost with low_k is {low_k_cost} while the cost ' +
    #            'with high k is {high_k_cost}').format(**vars()))

    subprocess.run(['dlclp'] + args + ['-K 0.001'])
    dlclp_info_file = open(fam_files_common_prefix + '.coal.klp.info', 'r')
    dlclp_runtime, low_dlclp_cost = get_runtime_and_cost(dlclp_info_file)
    dlclp_info_file.close()

    subprocess.run(['dlclp'] + args + ['-K 10'])
    dlclp_info_file = open(fam_files_common_prefix + '.coal.klp.info', 'r')
    dlclp_runtime, high_dlclp_cost = get_runtime_and_cost(dlclp_info_file)
    dlclp_info_file.close()

    # return low_k_cost, low_k_runtime, high_k_cost, high_k_runtime, dlclp_cost, dlclp_runtime
    return low_dlclp_cost, high_dlclp_cost


# =========
# main

# gather results from all families
with open(examples_dir + '/results4.tsv', 'w+') as log_file:
    # log_file.write('\t'.join(['sim', 'fam id', 'dlcpar and k = .001 cost', 'dlcpar and k = .001 runtime',
    #                           'dlcpar and k = 10 cost', 'dlcpar and k = 10 runtime', 'dlclp cost', 'dlclp runtime']) +
    #                '\n')

    log_file.write('\t'.join(['sim', 'fam id', 'dlclp low cost', 'dlclp high cost']) + '\n')

    # run each family
    for data_set, fam_id in sum([sorted_fams_in_data_set(data_set) for data_set in data_sets], []):
        print('fam: ', data_set, fam_id)
        # if data_set == '10e6-1x' or (data_set == '10e6-2x' and int(fam_id) < 231):
        #     continue

        prefix = '%s/%s/%s/%s' % (sim_flies_dir, data_set, fam_id, fam_id)
        result = test(prefix)

        log_file.write('\t'.join((data_set, fam_id) + result) + '\n')
        log_file.flush()
