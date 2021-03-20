from CleanLibraries import *
import matplotlib.pyplot as plt
import seaborn as sns

def get_basic_statistics_libraries(folder_path):
    depths = []
    read_lengths = []
    gcs = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith(".fastq"):
            path = folder_path + '/' + filename
            a,b,c = get_stats(path)
            
            depths.append(a)
            read_lengths.append(b)
            gcs.append(c)
            
    return depths, read_lengths, gcs
    
def plot_statistics(folder_path):
    depths, read_lengths, gcs = get_basic_statistics_libraries(folder_path)
    plt.figure(figsize=(15,10))
    plt.title('Read depth distribution among libraries', fontsize=20)
    sns.distplot(depths)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel(fontsize=20)
    plt.ylabel(fontsize=20)
    plt.show()
    
    plt.figure(figsize=(15,10))
    plt.title('Average read length distribution among libraries', fontsize=20)

    if np.unique(read_lengths).size==1:
        plt.hist(read_lengths, bins=len(read_lengths))
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(fontsize=20)
        plt.ylabel(fontsize=20)
        plt.show()
    else:
        sns.distplot(read_lengths)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(fontsize=20)
        plt.ylabel(fontsize=20)
        plt.show()

    plt.figure(figsize=(15,10))
    plt.title('Average %GC distribution among libraries', fontsize=20)
    plt.xlabel(fontsize=20)
    plt.ylabel(fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    sns.distplot(gcs)
    plt.show()
    

def get_values(lines):
    getvalues = False
    values = []

    for i, line in enumerate(lines):
        if(getvalues):
            vals = line.split('\t')
            values.append(float(vals[1][:-2]))
            if(int(vals[0]) == 100):
                return values

        if('normalized_position' in line):
            getvalues = True
            getvalues = True

def get_mapping_statistics(folder_path):
    df_positions = pd.DataFrame(columns=np.arange(0,101))
    count_table = pd.DataFrame()
    df_unmapped = pd.DataFrame(columns=['Unmapped'])
    df_read_depths = pd.DataFrame(columns=['ReadDepths'])
    df_mrna_bases = pd.DataFrame(columns=['mRNAbases'])
    
    j = 0

    for name in os.listdir(folder_path):
        if('.final.out' in name):
            w = open(folder_path + '/' + name)
            lines = w.readlines()
            unmapped_1 = float(lines[-4].split('|')[1][1:-2])
            unmapped_2 = float(lines[-6].split('|')[1][1:-2])
            unmapped_3 = float(lines[-8].split('|')[1][1:-2])
                              
            unmapped_percent = unmapped_1 + unmapped_2 + unmapped_3
            input_reads = int(lines[5].split('|')[1])
            
            df_read_depths.loc[name.split('.')[0]] = input_reads
            df_unmapped.loc[name.split('.')[0]] = unmapped_percent
            
            
        if('.RNA_Metrics' in name):
            w = open(folder_path + '/' + name)
            lines = np.array(w.readlines())
            df_positions.loc[name.split('.')[0]] = get_values(lines)
            if(len(lines) != 113):
            	print(len(lines))
            	continue
            df_mrna_bases.loc[name.split('.')[0]] = float(lines[7].split('\t')[21])
            j += 1
        
        if('fastqReadsPerGene' in name):
            gene_counts = pd.read_table(folder_path + '/' + name)
            gene_counts.columns = ['GeneID', 'Count1', 'Count2', 'Count3']
            gene_counts = gene_counts.iloc[4:, :2]
            gene_counts.columns = ['GeneID', name.split('.')[0]]
            
            try:
                count_table = count_table.merge(gene_counts, on='GeneID')
            except:
                count_table = gene_counts
                
    count_table.to_csv('CountTable.csv')
    return df_positions, count_table, df_read_depths, df_unmapped, df_mrna_bases


def plot_mapping_statistics(folder_path):
    df_pos, count_table, read_depths, unmapped_depths, mrna_bases = get_mapping_statistics(folder_path)
    plt.figure(figsize=(15,10))
    for i in range(len(df_pos)):
        plt.plot(df_pos.iloc[i])

    plt.title('Position Coverage', fontsize=20)
    plt.xlabel('Normalized Position', fontsize=20)
    plt.ylabel('Normalized Coverage',fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()
    
    genes = []
    input_reads = []
    unmapped = []
    
    for col in count_table:
        if(col != 'GeneID'):
            genes.append(count_table[col][count_table[col] != 0].shape[0])
            input_reads.append(read_depths.loc[col]['ReadDepths'])
            unmapped.append(unmapped_depths.loc[col]['Unmapped'])
           

    plt.figure(figsize=(15,10))
    g = sns.scatterplot(input_reads, genes)
    plt.title('Number of Genes per cell', fontsize=20)
    plt.xlabel('Number of input reads (million)', fontsize=20)
    plt.ylabel('Number of expressed genes', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    xlabels = [x for x in g.get_xticks()/1000000]
    g.set_xticklabels(xlabels)
    plt.show()
    
    plt.figure(figsize=(15,10))
    g = sns.scatterplot(input_reads, unmapped)
    plt.title('Fraction of unmapped reads per library', fontsize=20)
    plt.xlabel('Number of input reads (million)', fontsize=20)
    plt.ylabel('% Unmapped Reads', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    xlabels = [x for x in g.get_xticks()/1000000]
    g.set_xticklabels(xlabels)
    plt.show()
    
    plt.figure(figsize=(15,10))
    sns.distplot(mrna_bases['mRNAbases'], hist=False)
    plt.title('Fraction of mRNA reads per cell', fontsize=20)
    plt.xlabel('% mRNA bases', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()
    
    
if __name__ == "__main__":
    mode = sys.argv[1] 
    if(mode == 'basic'):
        plot_statistics(sys.argv[2])        
    elif(mode == 'mapping'):
        plot_mapping_statistics(sys.argv[2])
