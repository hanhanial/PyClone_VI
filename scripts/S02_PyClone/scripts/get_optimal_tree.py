import pandas as pd
import sys

def main(hf5file, sample, qip="QIP"):
    prefix = sample
    optimal_tree = pd.read_hdf(hf5file, key='results/optimal')
    optimal_tree = optimal_tree.loc[optimal_tree == True]
    for tree in optimal_tree.index:
        clonefreq = pd.read_hdf(
            hf5file, key='trees/{}/clone_freq'.format(tree))
        clonefreq['sample'] = clonefreq.index

        pd.melt(clonefreq, id_vars="sample", value_name="clonal_prev",
                var_name="clone_id").to_csv("{}_tree_{}_clone_freq.txt".format(prefix, tree),
                                            sep="\t", index=False)

        adjtree = pd.read_hdf(hf5file,
                              key='trees/{}/adjacency_list'.format(tree))
        adjtree.to_csv("{}_tree_{}.adj".format(prefix, tree), sep="\t", index=False,
                       header=False)

        if (qip=="QIP"):
            var_assignment = pd.read_hdf(
                hf5file, key='trees/{}/cluster_assignment'.format(tree))
            var_assignment.to_csv("{}_tree_{}_cluster_assignment.txt".format(prefix, tree), sep="\t",
                                index=False, header=False)
        else:
            var_assignment = pd.read_hdf(
                hf5file, key='trees/{}/variant_assignment'.format(tree))
            var_assignment.to_csv("{}_tree_{}_variant_assignment.txt".format(prefix, tree), sep="\t",
                                index=False, header=False)


if __name__ == "__main__":
    hf5file = sys.argv[1]
    sample = sys.argv[2]
    qip_ornot = sys.argv[3]
    main(hf5file, sample, qip=qip_ornot)
