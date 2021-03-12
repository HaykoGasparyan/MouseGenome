library(fastqcr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

fastqc_install()

fastqcFilesDir = args[1]
outputFilesDir = paste(args[1], '/FASTQCR', sep='')

fastqc(fq.dir = fastqcFilesDir,      # FASTQ files directory
       qc.dir = outputFilesDir ,     # Results directory
       threads = 25                
)

qc.dir = outputFilesDir
qc <- qc_aggregate(qc.dir)
statisticsTable = qc_stats(qc)
convertedTable = dcast(qc, paste("sample", "~", "module"), value.var = "status")

finalTable = merge(x = convertedTable, y = statisticsTable, by = "sample")

write.csv(finalTable, file = paste(outputFilesDir, '/Report.csv', sep=''))

