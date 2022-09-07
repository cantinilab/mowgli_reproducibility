# This file interprets the top peaks.

# Imports.
library(GenomicRanges)
library(motifmatchr)
library(chromVAR)
library(TFBSTools)
library(JASPAR2022)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(MuData)

# Read atac file.
in_atac <- "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/bottom_in_mofa.csv" # nolint
peaks_csv <- read.csv(in_atac, row.names = 2)

# Remove exotic chromosomes.
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "GL000194.1", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "GL000205.2", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "GL000205.2", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "GL000219.1", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "GL000219.1", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "KI270721.1", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "KI270726.1", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "KI270726.1", ]
peaks_csv <- peaks_csv[peaks_csv["Chromosome"] != "KI270713.1", ]

# Convert the peaks to GRanges.
chromosomes <- peaks_csv["Chromosome"][, 1]
ranges <- IRanges::IRanges(
    start = peaks_csv["Start"][, 1],
    end = peaks_csv["End"][, 1]
)
peaks <- GenomicRanges::GRanges(seqnames = chromosomes, ranges = ranges)

# Get JASPAR motifs.
opts <- list()
opts["species"] <- "Homo sapiens"
opts["collection"] <- "CORE"
motifs <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
motifs_pwm <- TFBSTools::toPWM(motifs)

# Get cisBP motifs.
data("human_pwms_v2")

# Fuse JASPAR and cisBP motifs.
for (name in names(motifs_pwm)) {
    human_pwms_v2[name] <- motifs_pwm[name]
}

# Create a Signac object from the peaks.
# Actually giving peaks_csv is nonsense.
# But we only care about the rownames so it's fine.
assay <- Signac::CreateChromatinAssay(
    peaks_csv,
    ranges = peaks,
    sep = c(":", "-")
)

# Create statistics about peaks.
assay <- Signac::RegionStats(
    object = assay,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# Add the downloaded motif PWM annotation.
assay <- Signac::AddMotifs(
    object = assay,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = human_pwms_v2
)

# Define where to save the motif enrichment outputs.
out_motif <- "/users/csb/huizing/Documents/PhD/Code/mowgli_reproducibility/enrich/bottom_motifs_mofa/motifs_" # nolint

# Get all bottom peaks.
background <- c()
for (dim in 0:14) {

    # Get the bottom peaks for that dimension.
    features <- rownames(assay)[peaks_csv[paste0("bottom_in_dim_", dim)] == "True"]

    background <- c(background, features)
}

# Iterate over MOFA+'s dimensions.
for (dim in 0:14) {

    # Get the bottom peaks for that dimension.
    features <- rownames(assay)[peaks_csv[paste0("bottom_in_dim_", dim)] == "True"]

    # Do motif enrichment analysis.
    enriched_motifs <- Signac::FindMotifs(
        object = assay,
        features = features,
        background = background #rownames(assay)
    )

    # Save the enrichment.
    write.csv(enriched_motifs, paste0(out_motif, dim, ".csv"))
}


