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
    plt.figure(figsize=(15,12))
    plt.title('Read depth distribution among libraries', fontsize=20)
    sns.distplot(depths)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()
    
    plt.figure(figsize=(15,12))
    plt.title('Average read length distribution among libraries', fontsize=20)

    if np.unique(read_lengths).size==1:
        plt.hist(read_lengths, bins=len(read_lengths))
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()
    else:
        sns.distplot(read_lengths)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()

    plt.figure(figsize=(15,12))
    plt.title('Average %GC distribution among libraries', fontsize=20)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    sns.distplot(gcs)
    plt.show()
    
    
if __name__ == "__main__":
    folder_path = sys.argv[1]
    plot_statistics(folder_path)





    