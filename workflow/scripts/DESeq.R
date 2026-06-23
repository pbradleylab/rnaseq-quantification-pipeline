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
			type="character", default="output.tsv"),
        make_option(c("-p", "--plot_path"),
                        type="character", default="output.png"),
        make_option(c("--plot_pdf_path"),
                        type="character", default="output.pdf"),
        make_option(c("--log2fc_threshold"),
                        type="double", default=0.6),
        make_option(c("--padj_threshold"),
                        type="double", default=0.05),
        make_option(c("--label_top_n"),
                        type="integer", default=10))

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)
dir.create(dirname(opt$output_file), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$plot_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$plot_pdf_path), recursive=TRUE, showWarnings=FALSE)

#Assign a variable to both files
count_data = read.csv(opt$counts_data, header = TRUE, sep = "\t")
metadata=read.csv(opt$metadata_file, header = TRUE)
fm = as.formula(paste0("~", opt$variable_to_analyze))

#We need to first change the first column into the row names with the following command
count_data = count_data %>% remove_rownames %>% column_to_rownames(var = "target_id")
if ("length" %in% colnames(count_data)) {
count_data = count_data[, colnames(count_data) != "length"]
}
#if it comes FALSE we need to turn it into a matrix using the following command, if it is TRUE the command is ommited
count_data_mtx = as.matrix(count_data)

#Very important to notice is that you can only have integer in your matrix. to ensure this the following command will round up all decimals

count_data_mtx = ceiling(count_data_mtx)

#Move the first col id names into the rownames so that later can match everything
sample_id_col = if ("sample_name" %in% colnames(metadata)) {
"sample_name"
} else if ("sample" %in% colnames(metadata)) {
"sample"
} else {
stop("Metadata must include a sample_name or sample column.")
}
metadata = metadata %>% remove_rownames %>% column_to_rownames(var = sample_id_col)


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
df$target_id = rownames(df)
df$padj_plot = ifelse(is.na(df$padj) | df$padj <= 0, NA, -log10(df$padj))
df$status = "Not significant"
df$status[
    !is.na(df$padj) &
        df$padj < opt$padj_threshold &
        df$log2FoldChange >= opt$log2fc_threshold
] = "Up"
df$status[
    !is.na(df$padj) &
        df$padj < opt$padj_threshold &
        df$log2FoldChange <= -opt$log2fc_threshold
] = "Down"

label_df = df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    arrange(padj) %>%
    head(opt$label_top_n)

plot = ggplot(df, aes(x=log2FoldChange, y=padj_plot, color=status)) +
    geom_point(alpha=0.75, na.rm=TRUE) +
    geom_text(
        data=label_df,
        aes(label=target_id),
        color="black",
        size=2.5,
        vjust=-0.6,
        check_overlap=TRUE,
        na.rm=TRUE
    ) +
    theme_minimal() +
    scale_color_manual(values=c("Down"="#2C7BB6", "Not significant"="#595959", "Up"="#D7191C")) +
    geom_vline(xintercept=c(-opt$log2fc_threshold, opt$log2fc_threshold), col="#D7191C", linetype="dashed") +
    geom_hline(yintercept=-log10(opt$padj_threshold), col="#D7191C", linetype="dashed") +
    labs(
        x="log2 fold change",
        y="-log10 adjusted p-value",
        color=NULL
    )

ggsave(opt$plot_path, plot=plot, width=7, height=5, dpi=300)
ggsave(opt$plot_pdf_path, plot=plot, width=7, height=5)
