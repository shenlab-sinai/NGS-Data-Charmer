import pandas as pd 
import sys

def generate_counts_matrix(outfile, input_files):
    counts_dict = {}
    for file in input_files:
        sample = file.split("/")[2].split(".counts.txt")[0]
        counts_dict[sample] = {}

        with open(file, "r") as infile:
            next(infile)
            next(infile)
            for lines in infile:
                lines = lines.strip().split("\t")
                counts_dict[sample][lines[0]] = int(float(lines[6]))

    dataframe = pd.DataFrame(counts_dict)
    dataframe.sort_index(inplace=True)
    dataframe.to_csv(outfile, sep='\t')

generate_counts_matrix(sys.argv[1], sys.argv[2:])
