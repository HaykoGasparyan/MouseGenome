import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import Counter
import re
import sys
import seaborn as sns
import os
from tqdm import tqdm_notebook as tqdm
from collections import defaultdict
import os
import numpy as np
import time


def parse_file(path):
    """
        Returns list of reads
    """
    
    
    with open(path, 'r') as f:
        content = f.readlines()
    return content


def list_duplicates(seq):
    """
        Returns each element and its frequency count for the given list
    """
    
    
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() if len(locs)>2)


def read_file(path):
    """
        Reads the fastq file for the given path
        Returns the list of original file lines and list of only read lines.
    """
    
    
    with open(path, 'r') as f:
        content = f.readlines()
        content_only_reads = np.array([line[:-1] for line in content[1::4]])
    return content, content_only_reads


def get_gc(read):
    """
        Returns GC% for the given read
    """
    
    GC = len("".join(re.split("[^C G]", read)))*100/len(read)
    return GC
    
    
def get_basic_stats(reads):
    lengths = []
    gcs =[]
    
    for line in reads:
        lengths.append(len(line))
        gcs.append(get_gc(line))
    
    return lengths, gcs
    
    
def get_n_count(read):
    """
        Returns N% for the read.
    """
    
    
    N_percent = len("".join(re.split("[^N]", read)))*100/len(read)
    return N_percent


def get_critical_values(file_reads):
    """
        Returns 0.25 and 0.75 quantiles of GCs and lengths of reads, and also N counts in reads
    """
    
    gcs = []
    ns = []
    lengths = []
    
    d_indexes = defaultdict(list)        
        
    for i,line in enumerate(file_reads):
        pass
        gcs.append(get_gc(line))
        ns.append(get_n_count(line))
        lengths.append(len(line))
        d_indexes[line].append(i)
        
    q1_gcs = np.percentile(gcs, 25)
    q3_gcs = np.percentile(gcs, 75)

    iqr_gcs = q3_gcs - q1_gcs
    
    g1 = q1_gcs - 1.5 * iqr_gcs
    g2 = q3_gcs + 1.5 * iqr_gcs
    
    q1_lengths = np.percentile(lengths, 25)
    q3_lengths = np.percentile(lengths, 75)

    iqr_lengths = q3_lengths - q1_lengths
    
    c1 = q1_lengths - 1.5 * iqr_lengths
    c2 = q3_lengths + 1.5 * iqr_lengths
    
    return g1, g2, gcs, c1, c2, lengths, ns, d_indexes


def remove_reads(file_original, file_reads):
    """
        Does basic preprocessing, removes bad reads
        ------------
        
        file_original: Original content of lines
        file_reads: Only read lines 
    """
    
    
    read_len_threshold = 30
    N_percent_threshold = 30
    
    new_file = []
    
    gc1, gc2, gcs, l1, l2, lengths, ns, dict_indexes = get_critical_values(file_reads)
    new_lengths = []
    new_gcs = []
    new_n_count = []

    for i in range(len(file_reads)):
        read = file_reads[i]
        indexes = dict_indexes[read][1:]
        
        GC = gcs[i]
        n_count = ns[i]
        length = lengths[i]
        
        if(length <= read_len_threshold):
            indexes.append(i)
        elif(length < l1):
            indexes.append(i)
        elif(GC < gc1 or GC > gc2):
            indexes.append(i)
        elif(n_count >= N_percent_threshold):
            indexes.append(i)
        else:
            new_gcs.append(GC)
            new_lengths.append(length)
            new_n_count.append(n_count)
        
        if(not i in indexes):
            new_file.append(file_original[4*i])
            new_file.append(file_original[4*i + 1])
            new_file.append(file_original[4*i + 2])
            new_file.append(file_original[4*i + 3])
            
    return new_file, new_gcs, new_lengths, new_n_count


def write_file(content, path):
    """
        Writes the given content into the given path and saves as .fastq file.
    """
    
    
    if not os.path.exists('/'.join(path.split('/')[:-1])):
        os.makedirs('/'.join(path.split('/')[:-1]))

    with open(path, 'w') as f:
        for line in content:
            f.write(line)
        f.close()
        

def preprocess_file(read_path, write_path):
    """
        Preproceessing of the given file.
        ------------
        read_path: File path to read
        write_path: Path to write the processed file
    """
    
    
    file_original, file_reads = read_file(read_path)
    new_file, gcs, lengths, n_counts = remove_reads(file_original, file_reads) 
        
    if(len(new_file)/4 <= 30_000):
        print('After processing Library depth is lower than 30_000. Removing Library {}...'.format(read_path.split('/')[-1]))
    elif(np.mean(lengths) <= 30):
        print('After processing Library average read lengths is lower than 30. Removing Library {}...'.format(read_path.split('/')[-1]))
    elif(np.mean(gcs) <= 30 or np.mean(gcs) >= 70):
        print('After processing Library average GC% is lower than 30 or higher than 70. Removing Library {}...'.format(read_path.split('/')[-1]))
    elif(np.mean(n_counts) >= 15):
        print('After processing Library average N% is higher than 15. Removing Library {}...'.format(read_path.split('/')[-1]))
    else:
        write_file(new_file, write_path)
        print('Cleared Library {}'.format(read_path.split('/')[-1]))
    
def get_stats(path):
    """
        Returns gc, depth and read lengths for the given libraries
    """
    
    file, file_original = read_file(path)
    lengths, gcs = get_basic_stats(file_original)
    depth = len(file_original)
    
    return depth, np.mean(lengths), np.mean(gcs)
    
def preprocess_libraries(folder_path):
    """
        Preprocessing of multiple libraries (database)
    """
    
    
    names = os.listdir(folder_path)
    for filename in tqdm(names):
        if('.fastq' in filename):
            _, file_reads = read_file(folder_path + '/' + filename)
            
            if(len(file_reads) <= 30_000):
                print('{} Library depth is lower than 30_000... Removing Library'.format(filename))
                continue
                
            preprocess_file(folder_path + '/' + filename, folder_path + '/' + 'Cleared' + '/' + filename)
            


if __name__ == "__main__":
    folder_path = sys.argv[1]
    preprocess_libraries(folder_path)










