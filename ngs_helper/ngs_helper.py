# HELPER FUNCTIONS

import os
import sys
from os import listdir
from os.path import isfile, join
import re
from Bio import SeqIO
import gzip
from collections import Counter

# Create the pattern of file endings

def create_endings(x):
    """
    Returns a list of likely fastq file endings

    Input Parameter: 
    x (int): Either 1 or 2, indicating the forward (1) or reverse (2) read.

    Returns: 
    list: A list of strings, representing the file endings a user might 
    use for denoting their fastq files.
    """
    return(["_R" + str(x) + "_001.fastq", "_R" + str(x) + "_001.fq",
            "_R" + str(x) + ".fastq", "_R" + str(x) + ".fq",
            "_" + str(x) + ".fastq", "_" + str(x) + ".fq",
            ".R" + str(x) + "_001.fastq", ".R" + str(x) + "_001.fq",
            ".R" + str(x) + ".fastq", ".R" + str(x) + ".fq",
            "." + str(x) + ".fastq", "." + str(x) + ".fq",
            "_r" + str(x) + "_001.fastq", "_r" + str(x) + "_001.fq",
            "_r" + str(x) + ".fastq", "_r" + str(x) + ".fq",
            ".r" + str(x) + "_001.fastq", ".r" + str(x) + "_001.fq",
            ".r" + str(x) + ".fastq", ".r" + str(x) + ".fq"])

# Function to list the fastq files present in the fastq folder


def getfilelist(myfastqpath):
    """
    Extracts fastq files from the files present in your fastq directory.

    Input Parameter: 
    myfastqpath (string): directory containing your fastq files.

    Returns:
    list: List containing two strings. 
    1st string is all non-metadata files in the fastq directory
    2nd string is all non-metadata files ending in '.gz'
    """
    onlyfiles = [f for f in listdir(myfastqpath) if
                 isfile(join(myfastqpath, f))]
    onlyfiles = [i for i in onlyfiles if
                 i.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))]
    gzfiles = [i for i in onlyfiles if i.endswith((".gz"))]
    return([onlyfiles, gzfiles])


def rename_files(oldname, replacement, myfastqpath):
    [os.rename(join(myfastqpath, i), join(myfastqpath, y))
     for i, y in zip(oldname, replacement)]

# Unify fastq files to single file ending


def fix_input_files(file_suffix, input_fileset, myfastqpath):
    """
    Renames mixed input fastq files to the most common file ending and 
    returns the selected file ending. NOTE: This step permenantly 
    renames your fastq files from their original file ending. 

    Input Parameter: 
    file_suffix (string): ".gz" or ""; Gzipped fastq files are expected 
    to end with the suffix ".gz". If files are NOT gzipped, the input is "".

    input_fileset (list): List of fastq file names to be examined. 
    As written, gzipped files are listed within the variable 'gzfiles' 
    and non-gzipped files are listed within the variable 'onlyfiles'.

    Returns: 
    list: A list containing four strings, the selected Read1 (forward read) 
    file ending and the corresponding Read2 (reverse read) file ending, 
    a list of all fastq-like files, and a list of gzipped fastq-like files.
    """

    # file_suffix, input_fileset = [".gz", gzfiles]
    # Create the series of fastq file endings to search
    base_endings_r1, base_endings_r2 = [create_endings(i) for i in (1, 2)]
    # Define the R1 and R2 suffix pairs for reference
    ending_dictionary = dict(zip(base_endings_r1, base_endings_r2))
    mylist = list()  # Create empty list

    # Traverse the R1 base endings to find the common ending
    for x in base_endings_r1:
        matched_ends = [
            i for i in input_fileset if i.endswith(x + file_suffix)]
        if(len(matched_ends) > 0):
            mylist.extend([x]*len(matched_ends))

    # If all samples are single-end
    if len(mylist) == 0:
        print("Your dataset appears to be entirely single-end files.")
        odd_files = [i for i in input_fileset
                     if i.endswith(".fq" + file_suffix)]
        if len(odd_files) > 0:
            old_rep = [i.replace(".fq" + suffix, ".fastq" + suffix)
                       for i in odd_files]
            rename_files(odd_files, old_rep, myfastqpath)

        # Re-assess fastq directory content and return filenames
        return([".fastq", ".fastq", getfilelist(myfastqpath)[0], 
            getfilelist(myfastqpath)[1]])

    # If R1 endings are present, check values and correct file names
    else:
        # create dictionary of mixed file endings
        mylist_endings = list(Counter(mylist).keys())
        # Find most common file ending
        myR1_suffix = max(Counter(mylist).items(), key=lambda x: x[1])[0]
        # Match chosen R1 ending to correct R2 ending
        myR2_suffix = ending_dictionary[myR1_suffix]
        # remove main R1 suffix from dictionary
        mylist_endings.remove(myR1_suffix)

        # Process forward reads
        if len(mylist_endings) > 0:
            for x in mylist_endings:
                oldnames = [
                    i for i in input_fileset if i.endswith(x + file_suffix)]
                old_rep = [i.replace(x, myR1_suffix) for i in oldnames]
                rename_files(oldnames, old_rep, myfastqpath)

            mylist = list()  # Create empty list to hold R2 file endings
            # Traverse the R2 base endings to endings
            for x in base_endings_r2:
                matched_ends = [i for i in input_fileset if i.endswith(
                    x + file_suffix) and x != myR2_suffix]
                if(len(matched_ends) > 0):
                    mylist.append(x)  # Create list of R2 files to be renamed
            if len(mylist) > 0: # Rename R2 files that don't match desired ending
                for x in mylist:
                    oldnames = [
                        i for i in input_fileset if i.endswith(x + file_suffix)]
                    old_rep = [i.replace(x, myR2_suffix) for i in oldnames]
                    rename_files(oldnames, old_rep, myfastqpath)

        # Re-assess file names
        if file_suffix == ".gz":
            input_fileset = getfilelist(myfastqpath)[1]
        else:
            input_fileset = getfilelist(myfastqpath)[0]

        # Now process single end files
        # Identify files that do not match the current R1, R2 ending
        odd_files = [i for i in input_fileset if not
                     i.endswith(myR1_suffix + file_suffix) if not
                     i.endswith(myR2_suffix + file_suffix)]

        # Partition single end files according to ending
        fastq_odd_1 = [i for i in odd_files if i.endswith(
            ".fastq" + file_suffix)]
        fastq_odd_2 = [i for i in odd_files if i.endswith(".fq" + file_suffix)]

        # If any apparently single-end files exist, then rename them
        if len(odd_files) > 0:
            print("Now unifying " + str(len(odd_files)) +
                  " single-end files to \"" + myR1_suffix +
                  file_suffix + "\" ending")
            # rename 'fastq' single-end files to correct ending
            if len(fastq_odd_1) > 0:
                old_rep = [i.replace(".fastq" + file_suffix,
                            myR1_suffix + file_suffix) for i in fastq_odd_1]
                rename_files(fastq_odd_1, old_rep, myfastqpath)
            # rename 'fq' single-end files to correct ending
            if len(fastq_odd_2) > 0:
                old_rep = [i.replace(".fq" + file_suffix,
                            myR1_suffix + file_suffix) for i in fastq_odd_2]
                rename_files(fastq_odd_2, old_rep, myfastqpath)
        # Re-assess and return filenames and file endings
        return([myR1_suffix, myR2_suffix, getfilelist(myfastqpath)[0], 
            getfilelist(myfastqpath)[1]])

# Function to retrieve and check cut&run read lengths


def check_readlength(suffix, input_fileset, R1_file_ending, myfastqpath):
    """
    When samples are specified to be cut&run:
    Samples the first read in each forward read file and 
    extracts the read length. 
    Checks if the samples have different read lengths.

    Input Parameter: 
    suffix (string): ".gz" or ""; Gzipped fastq files are expected to end 
    with the suffix ".gz". If files are NOT gzipped, the input is "".

    input_fileset (list): List of fastq file names to be examined. 
    As written, gzipped files are listed within the variable 'gzfiles' and 
    non-gzipped files are listed within the variable 'onlyfiles'.

    Returns:
    list: List containing two integers. 
    1st integer is the read length
    2nd integer is the read length, minus one
    """
    my_cr_files = [i for i in input_fileset if i.endswith(
        R1_file_ending + suffix)]
    if suffix == ".gz":
        try:  # Extract the read length of first forward read
            read_len_list = [len(next(SeqIO.parse(gzip.open(
                join(myfastqpath, i), "rt"), "fastq")).seq) for i in my_cr_files]
            dedup_lengths = list(dict.fromkeys(read_len_list))
        except:
            raise NameError(
                "One of your fastq files may be empty\nNow aborting...")
    elif suffix == "":
        try:
            read_len_list = [len(next(SeqIO.parse(join(myfastqpath, i),
                                "fastq")).seq) for i in my_cr_files]
            dedup_lengths = list(dict.fromkeys(read_len_list))
        except:
            raise NameError(
                "One of your fastq files may be empty\nNow aborting...")

    if len(dedup_lengths) > 1:
        raise NameError("Based on sampling the first read of each R1 fastq \
   file, your cut&run files have different read lengths!\nRecorded \
   lengths:" + ("elements in the list are " +
                ', '.join(['%.f']*len(dedup_lengths))) % tuple(dedup_lengths)
            + " base pairs" + "\nAborting...")
    else:
        print("Congratulations, your cut&run fastq files appear to have \
   uniform sequence lengths!\nProceeding with a read length of "
              + format(dedup_lengths[0]))
        return(dedup_lengths[0], dedup_lengths[0]-1)


# END HELPER FUNCTIONS
