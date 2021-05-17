import os
import sys
import pandas as pd

def main():

    # Set working directory to script's directory
    file_path = os.path.abspath(sys.argv[0])
    directory = os.path.dirname(file_path)
    os.chdir(directory)

    # feature table tsv is placed into a pandas dataframe
    fname = "16S_libraries_combined_taxonomy_feature_table"
    df = pd.read_csv(fname + ".tsv", sep="\t")

    # extract out all unassigned ASVs and place in a list
    unassigned_list = df[df["Taxon"] == "Unassigned"]["id"].tolist()

    # place all unassigned ASVs into a fasta file
    outfile = open(fname + "_unassigned.fasta", "w")
    for i in range(len(unassigned_list)):
        outfile.write(">unassigned_" + str(i) + "\n")
        outfile.write(unassigned_list[i] + "\n")
    outfile.close()

if __name__ == '__main__':
    main()