import numpy as np

def CHECK_CHR_SORTING(binIDs):
        
    chromosomes = [int(elem.split("_")[0]) for elem in binIDs.tolist()]
    bin_nums = [int(elem.split("_")[1]) for elem in binIDs.tolist()]
    
    unique_chromosomes = sorted([int(x) for x in list(set(chromosomes))])

    bool_chromosomes_list = list()
    bool_bin_nums_list = list()

    for chr in unique_chromosomes:

        idxs = np.array([idx for idx, elem in enumerate(chromosomes) if elem == chr])
        
        chromosomes_f = [chromosomes[idx] for idx in idxs]
        bin_nums_f = [bin_nums[idx] for idx in idxs]

        bool_chromosomes_list.append((chromosomes_f == sorted(chromosomes_f)))
        bool_bin_nums_list.append((bin_nums_f == sorted(bin_nums_f)))

    
    if (bool(np.all(bool_chromosomes_list)) and bool(np.all(bool_bin_nums_list))):
        return True
    else:
        return False