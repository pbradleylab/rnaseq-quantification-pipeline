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
        make_option(c("--plot_svg_path"),
                        type="character", default="output.svg"),
        make_option(c("--plot_pdf_path"),
                        type="character", default="output.pdf"),
        make_option(c("--expression_boxplot_path"),
                        type="character", default="normalized_expression_boxplot.png"),
        make_option(c("--expression_boxplot_svg_path"),
                        type="character", default="normalized_expression_boxplot.svg"),
        make_option(c("--expression_boxplot_pdf_path"),
                        type="character", default="normalized_expression_boxplot.pdf"),
        make_option(c("--expression_density_path"),
                        type="character", default="normalized_expression_density.png"),
        make_option(c("--expression_density_svg_path"),
                        type="character", default="normalized_expression_density.svg"),
        make_option(c("--expression_density_pdf_path"),
                        type="character", default="normalized_expression_density.pdf"),
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
dir.create(dirname(opt$plot_svg_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$plot_pdf_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$expression_boxplot_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$expression_boxplot_svg_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$expression_boxplot_pdf_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$expression_density_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$expression_density_svg_path), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(opt$expression_density_pdf_path), recursive=TRUE, showWarnings=FALSE)

#Assign a variable to both files
count_data = read.csv(opt$counts_data, header = TRUE, sep = "\t", check.names = FALSE)
metadata=read.csv(opt$metadata_file, header = TRUE, check.names = FALSE)
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


missing_counts = setdiff(rownames(metadata), colnames(count_data_mtx))
extra_counts = setdiff(colnames(count_data_mtx), rownames(metadata))
if (length(missing_counts) > 0 || length(extra_counts) > 0) {
stop(
    paste(
        "Count matrix columns and metadata samples do not match.",
        paste("Missing count columns:", paste(missing_counts, collapse=", ")),
        paste("Count columns without metadata:", paste(extra_counts, collapse=", ")),
        sep="\n"
    )
)
}
count_data_mtx = count_data_mtx[, rownames(metadata), drop=FALSE]

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
write.table(
    as.data.frame(res[order(res$padj),]),
    file=opt$output_file,
    sep="\t",
    quote=FALSE,
    col.names=NA
)
}

save_plot = function(plot, png_path, svg_path, pdf_path, width=7, height=5) {
ggsave(png_path, plot=plot, width=width, height=height, dpi=300)
svg(svg_path, width=width, height=height)
print(plot)
dev.off()
ggsave(pdf_path, plot=plot, width=width, height=height)
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

save_plot(plot, opt$plot_path, opt$plot_svg_path, opt$plot_pdf_path)

normalized_counts = counts(dds, normalized=TRUE)
normalized_df = as.data.frame(normalized_counts) %>%
    rownames_to_column(var="target_id") %>%
    pivot_longer(
        cols=-target_id,
        names_to="sample",
        values_to="normalized_count"
    ) %>%
    left_join(
        metadata %>% rownames_to_column(var="sample"),
        by="sample"
    ) %>%
    mutate(log10_normalized_count=log10(normalized_count + 1))

expression_boxplot = ggplot(
    normalized_df,
    aes(x=sample, y=log10_normalized_count, fill=.data[[opt$variable_to_analyze]])
) +
    geom_boxplot(outlier.size=0.4) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(
        x="sample",
        y="log10 normalized count + 1",
        fill=opt$variable_to_analyze
    )

expression_density = ggplot(
    normalized_df,
    aes(x=log10_normalized_count, color=.data[[opt$variable_to_analyze]])
) +
    geom_density(na.rm=TRUE) +
    theme_minimal() +
    labs(
        x="log10 normalized count + 1",
        y="density",
        color=opt$variable_to_analyze
    )

save_plot(
    expression_boxplot,
    opt$expression_boxplot_path,
    opt$expression_boxplot_svg_path,
    opt$expression_boxplot_pdf_path,
    width=8,
    height=5
)
save_plot(
    expression_density,
    opt$expression_density_path,
    opt$expression_density_svg_path,
    opt$expression_density_pdf_path
)
