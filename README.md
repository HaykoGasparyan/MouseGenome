# MouseGenome


1. Use ClearLibraries.py script to remove bad libraries and reads. It will create new directory ('Cleaned') in the same folder which will contain only needed things. You have to pass your folder (.fastq files are in) path as an argument


3. Use GetFastqcReport.r R script to run fastqc on your all libraries in the given folder and get Report.csv which will contain all aggregated FAILS and PASSES of your libraries by fastqc. You have to pass folder path as an argument.


4. Use PlotStatistics.py to get plots of basic statistics of your libraries. It will plot Depth distribution, average GC% content distribution and average read length distribution.You have to pass folder path as an argument.


Depending on your machine and size of data, this scripts may take some time, please be patient :)
