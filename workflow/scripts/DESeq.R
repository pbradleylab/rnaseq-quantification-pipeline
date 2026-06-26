library("DESeq2")
library("ggplot2")
library("tidyverse")
library("optparse")


option_list = list(
    make_option(c("--mode"), type="character", default="results"),
    make_option(c("-d", "--counts_data"), type="character", default=NULL),
    make_option(c("-m", "--metadata_file"), type="character", default=NULL),
    make_option(c("--design_formula"), type="character", default=NULL),
    make_option(c("-v", "--variable_to_analyze"), type="character", default=NULL),
    make_option(c("-r", "--reference_in_variable"), type="character", default=NULL),
    make_option(c("-o", "--output_file"), type="character", default=NULL),
    make_option(c("--svg_file"), type="character", default=NULL),
    make_option(c("--pdf_file"), type="character", default=NULL),
    make_option(c("--sample_report_file"), type="character", default=NULL),
    make_option(c("--transform_method"), type="character", default="vst"),
    make_option(c("--log2fc_threshold"), type="double", default=0.6),
    make_option(c("--padj_threshold"), type="double", default=0.05),
    make_option(c("--label_top_n"), type="integer", default=10)
)


parse_options = function() {
    opt_parser = OptionParser(option_list=option_list)
    opt = parse_args(opt_parser)
    if (is.null(opt$output_file)) {
        stop("--output_file is required.")
    }
    dir.create(dirname(opt$output_file), recursive=TRUE, showWarnings=FALSE)
    if (!is.null(opt$svg_file)) {
        dir.create(dirname(opt$svg_file), recursive=TRUE, showWarnings=FALSE)
    }
    if (!is.null(opt$pdf_file)) {
        dir.create(dirname(opt$pdf_file), recursive=TRUE, showWarnings=FALSE)
    }
    if (!is.null(opt$sample_report_file)) {
        dir.create(dirname(opt$sample_report_file), recursive=TRUE, showWarnings=FALSE)
    }
    opt
}


read_inputs = function(opt) {
    count_data = read.csv(opt$counts_data, header=TRUE, sep="\t", check.names=FALSE)
    metadata = read.csv(opt$metadata_file, header=TRUE, check.names=FALSE)

    count_data = count_data %>% remove_rownames %>% column_to_rownames(var="target_id")
    if ("length" %in% colnames(count_data)) {
        count_data = count_data[, colnames(count_data) != "length"]
    }
    count_data_mtx = ceiling(as.matrix(count_data))

    sample_id_col = if ("sample_name" %in% colnames(metadata)) {
        "sample_name"
    } else if ("sample" %in% colnames(metadata)) {
        "sample"
    } else {
        stop("Metadata must include a sample_name or sample column.")
    }
    metadata = metadata %>% remove_rownames %>% column_to_rownames(var=sample_id_col)

    if (!opt$variable_to_analyze %in% colnames(metadata)) {
        stop(paste("Metadata does not contain variable_to_analyze:", opt$variable_to_analyze))
    }
    if (!opt$reference_in_variable %in% metadata[[opt$variable_to_analyze]]) {
        stop(
            paste(
                "Metadata column",
                opt$variable_to_analyze,
                "does not contain reference level:",
                opt$reference_in_variable
            )
        )
    }
    metadata[[opt$variable_to_analyze]] = relevel(
        factor(metadata[[opt$variable_to_analyze]]),
        ref=opt$reference_in_variable
    )

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

    list(count_data_mtx=count_data_mtx, metadata=metadata)
}


fit_deseq = function(count_data_mtx, metadata, opt) {
    dds = DESeqDataSetFromMatrix(
        countData=count_data_mtx,
        colData=metadata,
        design=as.formula(opt$design_formula)
    )
    DESeq(dds)
}


get_contrast_results = function(dds, opt) {
    contrast_levels = setdiff(
        levels(colData(dds)[[opt$variable_to_analyze]]),
        opt$reference_in_variable
    )
    if (length(contrast_levels) < 1) {
        stop(
            paste(
                "Metadata column",
                opt$variable_to_analyze,
                "must contain at least one non-reference level."
            )
        )
    }
    contrast_level = contrast_levels[1]
    results(dds, contrast=c(opt$variable_to_analyze, contrast_level, opt$reference_in_variable))
}


format_results = function(res, opt) {
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
    df
}


save_table = function(res, opt) {
    write.table(
        as.data.frame(res[order(res$padj),]),
        file=opt$output_file,
        sep="\t",
        quote=FALSE,
        col.names=NA
    )
}


save_matrix = function(matrix_data, output_file) {
    df = as.data.frame(matrix_data, check.names=FALSE) %>%
        rownames_to_column(var="target_id")
    write.table(
        df,
        file=output_file,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
}


save_plot = function(plot, opt, width=7, height=5) {
    ggsave(opt$output_file, plot=plot, width=width, height=height, dpi=300)
    if (!is.null(opt$svg_file)) {
        svg(opt$svg_file, width=width, height=height)
        print(plot)
        dev.off()
    }
    if (!is.null(opt$pdf_file)) {
        ggsave(opt$pdf_file, plot=plot, width=width, height=height)
    }
}


get_transformed_dds = function(dds, opt) {
    method = tolower(opt$transform_method)
    if (method == "rlog") {
        return(rlog(dds, blind=FALSE))
    }
    if (method == "auto") {
        if (ncol(dds) > 30) {
            method = "vst"
        } else {
            return(rlog(dds, blind=FALSE))
        }
    }
    if (method != "vst") {
        stop(paste("Unsupported transform_method:", opt$transform_method))
    }
    if (nrow(dds) < 1000) {
        varianceStabilizingTransformation(dds, blind=FALSE)
    } else {
        vst(dds, blind=FALSE)
    }
}


make_cooks_reports = function(dds, res, opt) {
    cooks = assays(dds)[["cooks"]]
    if (is.null(cooks)) {
        stop("DESeq2 object does not contain a Cook's distance assay.")
    }
    model_matrix = model.matrix(design(dds), colData(dds))
    cooks_cutoff = NA_real_
    residual_df = ncol(dds) - ncol(model_matrix)
    if (residual_df > 0) {
        cooks_cutoff = qf(0.99, ncol(model_matrix), residual_df)
    }
    max_cooks = apply(cooks, 1, max, na.rm=TRUE)
    mean_cooks = rowMeans(cooks, na.rm=TRUE)
    cooks_outlier = if (is.na(cooks_cutoff)) {
        rep(NA, length(max_cooks))
    } else {
        max_cooks > cooks_cutoff
    }
    sample_with_max_cooks = colnames(cooks)[max.col(cooks, ties.method="first")]
    result_df = as.data.frame(res)
    gene_report = data.frame(
        target_id=rownames(cooks),
        baseMean=result_df[rownames(cooks), "baseMean"],
        pvalue=result_df[rownames(cooks), "pvalue"],
        padj=result_df[rownames(cooks), "padj"],
        max_cooks=max_cooks,
        mean_cooks=mean_cooks,
        sample_with_max_cooks=sample_with_max_cooks,
        cooks_cutoff=cooks_cutoff,
        cooks_outlier=cooks_outlier,
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    sample_report = data.frame(
        sample=colnames(cooks),
        max_cooks=apply(cooks, 2, max, na.rm=TRUE),
        mean_cooks=colMeans(cooks, na.rm=TRUE),
        genes_over_cooks_cutoff=if (is.na(cooks_cutoff)) {
            NA
        } else {
            colSums(cooks > cooks_cutoff, na.rm=TRUE)
        },
        cooks_cutoff=cooks_cutoff,
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    write.table(
        gene_report[order(gene_report$max_cooks, decreasing=TRUE),],
        file=opt$output_file,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
    if (!is.null(opt$sample_report_file)) {
        write.table(
            sample_report[order(sample_report$max_cooks, decreasing=TRUE),],
            file=opt$sample_report_file,
            sep="\t",
            quote=FALSE,
            row.names=FALSE
        )
    }
}


make_volcano_plot = function(df, opt) {
    label_df = df %>%
        filter(!is.na(padj), !is.na(log2FoldChange)) %>%
        arrange(padj) %>%
        head(opt$label_top_n)

    ggplot(df, aes(x=log2FoldChange, y=padj_plot, color=status)) +
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
        labs(x="log2 fold change", y="-log10 adjusted p-value", color=NULL)
}


make_ma_plot = function(df, opt) {
    ggplot(df, aes(x=baseMean, y=log2FoldChange, color=status)) +
        geom_point(alpha=0.75, na.rm=TRUE) +
        scale_x_log10() +
        scale_color_manual(values=c("Down"="#2C7BB6", "Not significant"="#595959", "Up"="#D7191C")) +
        geom_hline(yintercept=c(-opt$log2fc_threshold, opt$log2fc_threshold), col="#D7191C", linetype="dashed") +
        geom_hline(yintercept=0, col="#595959", linetype="dotted") +
        theme_minimal() +
        labs(x="mean normalized expression", y="log2 fold change", color=NULL)
}


make_expression_boxplot = function(dds, metadata, opt) {
    normalized_df = as.data.frame(counts(dds, normalized=TRUE)) %>%
        rownames_to_column(var="target_id") %>%
        pivot_longer(cols=-target_id, names_to="sample", values_to="normalized_count") %>%
        left_join(metadata %>% rownames_to_column(var="sample"), by="sample") %>%
        mutate(log10_normalized_count=log10(normalized_count + 1))

    ggplot(
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
}


make_expression_density = function(dds, metadata, opt) {
    normalized_df = as.data.frame(counts(dds, normalized=TRUE)) %>%
        rownames_to_column(var="target_id") %>%
        pivot_longer(cols=-target_id, names_to="sample", values_to="normalized_count") %>%
        left_join(metadata %>% rownames_to_column(var="sample"), by="sample") %>%
        mutate(log10_normalized_count=log10(normalized_count + 1))

    ggplot(
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
}


make_sample_distance_heatmap = function(vsd) {
    sample_distances = as.matrix(dist(t(assay(vsd))))
    sample_distance_df = as.data.frame(as.table(sample_distances))
    colnames(sample_distance_df) = c("sample_a", "sample_b", "distance")
    sample_distance_df$sample_a = factor(
        sample_distance_df$sample_a,
        levels=colnames(sample_distances)
    )
    sample_distance_df$sample_b = factor(
        sample_distance_df$sample_b,
        levels=rev(colnames(sample_distances))
    )

    ggplot(sample_distance_df, aes(x=sample_a, y=sample_b, fill=distance)) +
        geom_tile(color="white", linewidth=0.2) +
        coord_fixed() +
        scale_fill_gradient(low="#F7FBFF", high="#08306B") +
        theme_minimal() +
        theme(
            axis.text.x=element_text(angle=45, hjust=1),
            panel.grid=element_blank()
        ) +
        labs(x=NULL, y=NULL, fill="distance")
}


make_pca_plot = function(vsd, metadata, opt) {
    pca = prcomp(t(assay(vsd)))
    percent_var = pca$sdev^2 / sum(pca$sdev^2)
    pca_df = as.data.frame(pca$x[, c("PC1", "PC2"), drop=FALSE]) %>%
        rownames_to_column(var="sample") %>%
        left_join(metadata %>% rownames_to_column(var="sample"), by="sample")

    ggplot(
        pca_df,
        aes(x=PC1, y=PC2, color=.data[[opt$variable_to_analyze]], label=sample)
    ) +
        geom_point(size=3, alpha=0.85) +
        geom_text(color="black", size=3, vjust=-0.7, check_overlap=TRUE) +
        theme_minimal() +
        labs(
            x=paste0("PC1: ", round(percent_var[1] * 100), "% variance"),
            y=paste0("PC2: ", round(percent_var[2] * 100), "% variance"),
            color=opt$variable_to_analyze
        )
}


make_library_size_factor_plot = function(count_data_mtx, dds, metadata, opt) {
    library_qc_df = data.frame(
        sample=colnames(count_data_mtx),
        raw_library_size=colSums(count_data_mtx),
        deseq2_size_factor=sizeFactors(dds),
        check.names=FALSE
    ) %>%
        left_join(metadata %>% rownames_to_column(var="sample"), by="sample") %>%
        pivot_longer(
            cols=c(raw_library_size, deseq2_size_factor),
            names_to="metric",
            values_to="value"
        ) %>%
        mutate(
            metric=recode(
                metric,
                raw_library_size="raw library size",
                deseq2_size_factor="DESeq2 size factor"
            )
        )

    ggplot(
        library_qc_df,
        aes(x=sample, y=value, fill=.data[[opt$variable_to_analyze]])
    ) +
        geom_col(width=0.75) +
        facet_wrap(~metric, scales="free_y", ncol=1) +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        labs(x="sample", y=NULL, fill=opt$variable_to_analyze)
}


run = function() {
    opt = parse_options()
    inputs = read_inputs(opt)
    dds = fit_deseq(inputs$count_data_mtx, inputs$metadata, opt)
    res = get_contrast_results(dds, opt)

    if (opt$mode == "results") {
        save_table(res, opt)
    } else if (opt$mode == "normalized_counts") {
        save_matrix(counts(dds, normalized=TRUE), opt$output_file)
    } else if (opt$mode == "transformed_counts") {
        save_matrix(assay(get_transformed_dds(dds, opt)), opt$output_file)
    } else if (opt$mode == "cooks_report") {
        make_cooks_reports(dds, res, opt)
    } else if (opt$mode == "volcano") {
        save_plot(make_volcano_plot(format_results(res, opt), opt), opt)
    } else if (opt$mode == "ma") {
        save_plot(make_ma_plot(format_results(res, opt), opt), opt)
    } else if (opt$mode == "expression_boxplot") {
        save_plot(make_expression_boxplot(dds, inputs$metadata, opt), opt, width=8, height=5)
    } else if (opt$mode == "expression_density") {
        save_plot(make_expression_density(dds, inputs$metadata, opt), opt)
    } else if (opt$mode == "sample_distance_heatmap") {
        save_plot(make_sample_distance_heatmap(get_transformed_dds(dds, opt)), opt, width=7, height=6)
    } else if (opt$mode == "pca") {
        save_plot(make_pca_plot(get_transformed_dds(dds, opt), inputs$metadata, opt), opt)
    } else if (opt$mode == "library_size_factors") {
        save_plot(
            make_library_size_factor_plot(inputs$count_data_mtx, dds, inputs$metadata, opt),
            opt,
            width=8,
            height=6
        )
    } else {
        stop(paste("Unsupported DESeq mode:", opt$mode))
    }
}


run()
