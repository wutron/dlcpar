
import os, sys, glob
import subprocess


def test_dlclp_generates_correct_recon(stree_path, smap_path, fam_files_common_prefix):
    subprocess.run(['dlclp', '-s', stree_path, '-S', smap_path, '-I', '.coal.tree',
                    fam_files_common_prefix + '.coal.tree'])
    subprocess.run(['dlcpar_to_dlcoal', '-s', stree_path, fam_files_common_prefix + '.dlclp.tree'])
    command = ['dlcoal_eq', '-s', stree_path, '-S', smap_path, '--use-locus-mpr',
               fam_files_common_prefix, fam_files_common_prefix + '.dlclp']
    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    process.wait()
    out, err = process.communicate()
    print(out)
    print('error: ' + str(err))

    # cleanup: delete all the files that were outputed by these commands
    for output_file in glob.glob(fam_files_common_prefix + '.dlclp*'):
        os.remove(output_file)

    return 'True' in str(out)


def main():
    with open('/research/yjw/jgardi/work/dlcpar_ilp/config/test_fly_fams.tsv', 'r') as test_fams:
        index_in_test_fams = int(sys.argv[1])
        line_for_fam = test_fams.read().splitlines()[index_in_test_fams]
        stree_path, smap_path, dataset_path, fam_id = line_for_fam.split('\t')

        fam_files_common_prefix = '%s/%s/%s' % (dataset_path, fam_id, fam_id)

        if test_dlclp_generates_correct_recon(stree_path, smap_path, fam_files_common_prefix):
            print(dataset_path + '/' + fam_id, 'succeeded')
        else:
            print(dataset_path + '/' + fam_id, 'failed')


if __name__ == '__main__':
    main()
