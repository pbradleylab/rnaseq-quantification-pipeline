#First open the programs to use
library("DESeq2")
library("ggplot2")
library("tidyverse")
library("optparse")
if (sys.nframe() == 0) {
option_list = list(
	make_option(c("-d", "--counts_data"),
			type="character", default=NULL),
	make_option(c("-m", "--metadata_file"),
			type="character", default=NULL),
	make_option(c("-v", "--variable_to_analyze"),
			type="character", default=NULL),
	make_option(c("-r", "--reference_in_variable"),
			type="character", default=NULL),
	make_option(c("-o", "--output_file"),
			type="character", default="output.tsv"))

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

#Assign a variable to both files
count_data = read.csv(opt$counts_data, header = TRUE)
metadata=read.csv(opt$metadata_file, header = TRUE)
fm = as.formula(paste0("~", opt$variable_to_analyze))

#We need to first change the first column into the row names with the following command
count_data = count_data %>% remove_rownames %>% column_to_rownames(var = "target_id")
count_data = count_data[,-c(1)]
#if it comes FALSE we need to turn it into a matrix using the following command, if it is TRUE the command is ommited 
count_data_mtx = as.matrix(count_data)

#Very important to notice is that you can only have integer in your matrix. to ensure this the following command will round up all decimals

count_data_mtx = ceiling(count_data_mtx)

#Move the first col id names into the rownames so that later can match everything
metadata = metadata %>% remove_rownames %>% column_to_rownames(var= "sample_name")


#Now we verify if the names from the data and metadata matches (both must come TRUE)
all(rownames(metadata) %in% colnames(count_data_mtx))

all(rownames(metadata) == colnames(count_data_mtx))

#now we can run the DESeq with our matrix MAKE SURE TO CHANGE THE DESIGN TO THE NAME OF THE SECOND COLUMN FROM THE METADATA)
dds = DESeqDataSetFromMatrix(countData = count_data_mtx, colData = metadata, design = fm)

#now we need to set the reference for our study MAKE SURE YOU CHANGE THE ref TO ONE OF YOUR CONDITIONS IN THE COLUMN TWO OF METADATA

dds[[opt$variable_to_analyze]] = relevel(dds[[opt$variable_to_analyze]], ref = opt$reference_in_variable)

#now we can run the differential gene expression analysis

dds = DESeq(dds)

#Still trying to find what the following command does but it need to be done
res = results(dds)

#now we need to order the gene expression by adjusted p-value
res[order(res$padj),]

#finally you can import the data
write.csv(as.data.frame(res[order(res$padj),]), file = opt$output_file)
}

#Read and assign a variable to the output file
df = as.data.frame(res[order(res$padj),])

#Generate a preliminary graph for the user to see
plot = ggplot(df, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.05), col="red")
png("volcano_plot.png")
print(plot)
dev.off()
