library(data.table)
library(ggplot2)
library(qqman)

gwasdir = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/gwas_results/"
plotdir = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/plots/"
PHENOS = c("gfap", "nfl", "ptau", "totaltau", "abeta_ratio")

for(i in 1:5){
    pheno = PHENOS[i]
    gwas <- fread(
        file.path(gwasdir, paste0(pheno, ".assoc.txt")),
        data.table = FALSE
    )

    gwas$chr <- as.integer(gwas$chr)
    gwas$ps  <- as.integer(gwas$ps)
    gwas$p_wald <- as.numeric(gwas$p_wald)

    gwas <- gwas[!is.na(gwas$p_wald) & gwas$p_wald >= 0 & gwas$p_wald <= 1, ]

    n <- nrow(gwas)
    qq_df <- data.frame(
        exp = -log10(ppoints(n)),
        obs = -log10(sort(gwas$p_wald))
    )
    chisq <- qchisq(1 - gwas$p_wald, df = 1)
    lambda <- median(chisq) / qchisq(0.5, 1)

    ggplot(qq_df, aes(x = exp, y = obs)) +
        geom_point(size = 0.6, alpha = 0.6) +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        labs(
            title = paste0(pheno, " QQ plot (λ = ", round(lambda, 3), ")"),
            x = expression(Expected~~-log[10](p)),
            y = expression(Observed~~-log[10](p))
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank()
        )
    ggsave(
        filename = file.path(plotdir, paste0(pheno, "_QQplot.png")),
        width = 6,
        height = 6,
        dpi = 720
    )

    manh <- data.frame(
        CHR = gwas$chr,
        BP  = gwas$ps,
        P   = gwas$p_wald,
        SNP = gwas$rs
    )

    png(file.path(plotdir, paste0(pheno, "_Manhattan.png")), width = 1600, height = 800)
    manhattan(
        manh,
        main = paste(pheno, "Manhattan plot"),
        genomewideline = -log10(5e-8),
        suggestiveline = -log10(1e-5)
    )
    dev.off()
}