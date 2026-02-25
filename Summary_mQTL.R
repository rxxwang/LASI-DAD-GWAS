library(data.table)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20260113_GWAS1/")
result <- fread("mqtl_summary_report.txt")

library(ggplot2)
library(scales) # for scientific notation

library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# --- 1. Settings ---
cpg_id <- "cg05615147"  # Replace with your actual CpG name
chr_num <- 15           # Replace with actual chromosome
cpg_pos <- 83067164     # Replace with actual position
window <- 1000000
results_path <- paste0("mQTL_results/", cpg_id, ".assoc.txt")

# --- 2. Load Data ---
# Column 2 = rsID, Column 3 = pos, Column 12 = p_wald
dt <- fread(results_path)

# 2. Create the Data Track (The Manhattan Plot)
# We convert your GEMMA results into a GRanges object
chr_label = paste0("chr", chr_num)
dTrack <- DataTrack(range = GRanges(seqnames = chr_label, 
                                    ranges = IRanges(start = dt$ps, end = dt$ps)),
                    data = -log10(dt$p_wald), 
                    name = "-log10(P)",
                    type = "p", # 'p' for points
                    col = "royalblue",
                    cex = 0.8)

# 3. Create the Gene Track (Exons/Introns)
# This pulls automatically from the UCSC database
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gTrack <- GeneRegionTrack(txdb, chromosome = chr_label, 
                          start = cpg_pos - window, 
                          end = cpg_pos + window,
                          showId = TRUE, # Shows gene names
                          geneSymbol = TRUE,
                          name = "Genes",
                          transcriptAnnotation = "symbol",
                          background.title = "brown")

# 4. Create an Axis Track
aTrack <- GenomeAxisTrack()

# 5. Highlight the CpG site
ht <- HighlightTrack(trackList = list(dTrack, gTrack),
                     start = cpg_pos, end = cpg_pos,
                     chromosome = chr_label,
                     col = "red", fill = "#FFE6E6")

# 6. Plot everything together
png(filename = paste0(cpg_id, "_genomic_region.png"), 
    width = 10, height = 7, units = "in", res = 300)
plotTracks(list(aTrack, ht), 
           from = cpg_pos - window, 
           to = cpg_pos + window,
           main = paste("Regional mQTL & Gene Map:", cpg_id),
           cex.main = 1.2)
dev.off()


library(org.Hs.eg.db) # For gene type mapping

# --- 1. Filter for Protein Coding Genes ---
# Get all Entrez IDs that are protein coding
all_genes <- keys(org.Hs.eg.db, keytype = "GENETYPE")
protein_coding_ids <- select(org.Hs.eg.db, 
                             keys = "protein-coding", 
                             keytype = "GENETYPE", 
                             columns = "ENTREZID")$ENTREZID

# --- 2. Create the Gene Track ---
# We use the 'geneSymbol = TRUE' to map Entrez IDs to names like "TP53"
gTrack <- GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, 
                          chromosome = chr_label, 
                          start = cpg_pos - window, 
                          end = cpg_pos + window,
                          name = "Protein Coding Genes",
                          transcriptAnnotation = "symbol", # Use Gene Name
                          collapseTranscripts = "meta",     # Merges isoforms into 1 line
                          fill = "darkorange",
                          col = "black")

# --- 3. Filter the track data specifically ---
# This removes non-protein-coding genes from the track before plotting
feature(gTrack) <- "protein_coding" # Cosmetic label
ranges(gTrack) <- ranges(gTrack)[gene(gTrack) %in% protein_coding_ids]

# --- 4. Plot and Save ---
png(filename = paste0(cpg_id, "_protein_coding_region.png"), 
    width = 10, height = 6, units = "in", res = 300)

plotTracks(list(aTrack, ht), 
           from = cpg_pos - window, 
           to = cpg_pos + window,
           main = paste("mQTL Regional Map (Protein Coding):", cpg_id),
           cex.main = 1,
           # The following settings clean up the labels
           showExonId = FALSE,      # Remove exon names
           just.group = "right",    # Put gene name to the right of the gene
           geneSymbols = TRUE,      # Ensure symbols are used
           fontcolor.group = "black", 
           fontsize.group = 10)

dev.off()