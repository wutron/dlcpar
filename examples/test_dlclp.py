import re
import subprocess
import sys
import os
import glob


def get_runtime_and_cost(file):
    lines = file.readlines()
    runtime_line = [x for x in lines if 'Runtime:' in x][0]
    cost_line = [x for x in lines if 'Cost:' in x][0]

    return re.search(r'\d+\.\d+', runtime_line).group(), re.search(r'\d+\.\d+', cost_line).group()


def test(stree_path, smap_path, fam_files_common_prefix):
    """
    :param str fam_files_common_prefix: the prefix sahred by all the files for this family
    :return tuple: the costs and runtime from testing on this family
    """
    coal_tree = fam_files_common_prefix + '.coal.tree'
    args = ['-s', stree_path, '-S', smap_path, '-I',  '.coal.tree', coal_tree]
    dlcpar_base_command = ['dlcpar'] + args

    subprocess.run(dlcpar_base_command + ['-K 0.001'])
    dlcpar_info_file = open(fam_files_common_prefix + '.dlcpar.info', 'r')
    low_k_runtime, low_k_cost = get_runtime_and_cost(dlcpar_info_file)
    dlcpar_info_file.close()

    subprocess.run(dlcpar_base_command + ['-K 10'])
    dlcpar_info_file = open(fam_files_common_prefix + '.dlcpar.info', 'r')
    high_k_runtime, high_k_cost = get_runtime_and_cost(dlcpar_info_file)
    dlcpar_info_file.close()

    dlclp_base_command = ['dlclp'] + args

    subprocess.run(dlclp_base_command + ['-K 0.001'])
    dlclp_info_file = open(fam_files_common_prefix + '.dlclp.info', 'r')
    dlclp_low_k_runtime, dlclp_low_k_cost = get_runtime_and_cost(dlclp_info_file)
    dlclp_info_file.close()

    subprocess.run(dlclp_base_command + ['-K 10'])
    dlclp_info_file = open(fam_files_common_prefix + '.dlclp.info', 'r')
    dlclp_high_k_runtime, dlclp_high_k_cost = get_runtime_and_cost(dlclp_info_file)
    dlclp_info_file.close()

    # cleanup: delete all the files that were outputed by these commands
    for output_file in glob.glob(fam_files_common_prefix + '.dlclp*') + glob.glob(fam_files_common_prefix + '.dlcpar*'):
        os.remove(output_file)

    return low_k_cost, low_k_runtime, high_k_cost, high_k_runtime, dlclp_low_k_cost, dlclp_low_k_runtime, \
           dlclp_high_k_cost, dlclp_high_k_runtime


def main():
    results_path = '/research/yjw/jgardi/work/dlcpar_ilp/analysis/results.tsv'
    test_fam_path = '/research/yjw/jgardi/work/dlcpar_ilp/config/test_fly_fams.tsv'
    with open(results_path, 'a+') as log_file, open(test_fam_path, 'r') as test_fams:
        index_in_test_fams = int(sys.argv[1])
        line_for_fam = test_fams.read().splitlines()[index_in_test_fams]
        stree_path, smap_path, dataset_path, fam_id = line_for_fam.split('\t')

        fam_files_common_prefix = '%s/%s/%s' % (dataset_path, fam_id, fam_id)
        result = test(stree_path, smap_path, fam_files_common_prefix)

        log_file.write('\t'.join((dataset_path, fam_id) + result) + '\n')


if __name__ == '__main__':
    main()
