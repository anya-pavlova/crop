import numpy as np
import pandas as pd
import sys
import getopt
import os

ppm = 10
pos_add = {'M+H': 1.007276, 'M+NH4': 18.033823, 'M+Na': 22.989218}


def usage():
    print("")

def process_groups(groups_file_name, database_df, output_file_name):
    groups_df = pd.read_csv(groups_file_name, sep="\t",
                             header=0, index_col=0, low_memory=False)
    
    
    masses = database_df["exact_mass"]
    labels = database_df["lm_id"].values  
    indices = groups_df.index 

    annotations_indices = []
    annotations_labels = []
    annotations_adducts = [] 

    for adduct, adduct_mass in pos_add.items():
        # Idea: |mz-m|/m < delta <---> mz/(1+delta) < m < mz/(1-delta)
        # Use binary search to accelerate search

        mz = groups_df['mzmed']-adduct_mass
            
        delta = ppm/1000000
        
        bottom_bounds = mz/(1+delta)
        upper_bounds = mz/(1-delta)
        
        bottom_indices = masses.searchsorted(bottom_bounds, side='right')
        upper_indices = masses.searchsorted(upper_bounds, side='left')

        for index, bottom_ind, upper_ind in zip(indices, bottom_indices, upper_indices):
            k = upper_ind-bottom_ind

            if k==0:
                annotations_indices.append(index)
                annotations_labels.append(np.nan)
                annotations_adducts.append(adduct)
            else:
                annotations_indices.extend([index]*k)
                annotations_labels.extend(labels[bottom_ind:upper_ind])
                annotations_adducts.extend([adduct]*k)

    annotations = pd.DataFrame({'lm_id':annotations_labels, 'adduct':annotations_adducts}, index=annotations_indices)
    
    os.makedirs(os.path.dirname(output_file_name), exist_ok = True)
    annotations.to_csv(output_file_name)


def main(argv):
    annotation_db = ''
    input_directory = ''
    output_directory = ''

    try:
        opts, args = getopt.getopt(argv,"hd:i:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-d"):
            annotation_db = arg
        elif opt in ("-i"):
            input_directory = arg
        elif opt in ("-o"):
            output_directory = arg

    database_df = pd.read_csv(annotation_db, sep="\t",
                           header=0,
                           index_col=None,
                           low_memory=False)
    
    database_df.sort_values("exact_mass", inplace=True)

    for filename in os.listdir(input_directory):
        process_groups(os.path.join(input_directory, filename), database_df, os.path.join(output_directory, filename + ".ann.txt"))

if __name__ == "__main__":
    main(sys.argv[1:])
