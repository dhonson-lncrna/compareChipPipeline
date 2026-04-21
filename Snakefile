'''
A Snakemake pipeline for batch aware differential peak calling
'''

import os
from itertools import combinations

##############################################################################
# Load required variables
##############################################################################

configfile: './config.yaml'

# Get input directories
BAM_DIR = config.get('bam_dir')
MACS3_DIR = config.get('macs3_dir')
LFC_BW = config.get('lfc_bw_dir')

# Make output directory
DIR_OUT = config.get('dir_out')
os.makedirs(DIR_OUT, exist_ok=True)
    
SUBDIRS = ['counts',
           'countMatrix',
           'deseq2',
           'mergePeaks',
           'avg_lfc',
           'intersect',
           'heatmap']
DIR_DICT = {d: os.path.join(DIR_OUT, d) for d in SUBDIRS}

# Parameters
BAM_SUFFIX = config.get('bam_suffix')
MIN_OVERLAP = config.get('min_overlap')
ALPHA = config.get('alpha')
LFC_THRESHOLD = config.get('lfc_threshold')
GENE_BODIES_BED = config.get('gene_bodies_bed')
DT_A = config.get('a')
DT_B = config.get('b')
DT_M = config.get('m')
DT_BS = config.get('bs')

##############################################################################
# Generate BAM file metadata
##############################################################################

BAM_FILES = [f for f in os.listdir(BAM_DIR) if f.endswith(BAM_SUFFIX)]
BAM_FILES = [f for f in BAM_FILES if 'input' not in f]
MARKS = set([i.split('_')[1] for i in BAM_FILES])

REPS = set([i.split('-')[0] for i in BAM_FILES])
REPJOIN = '-'.join(sorted(REPS))

#CONDITIONS = set([i.split('-')[1].split('_')[0] for i in BAM_FILES])
CONDITIONS = ["apo", "wt"]
CONDITIONS_LO = [f"{c}-lo" for c in CONDITIONS]

BIGWIGS = expand(
    os.path.join(DIR_DICT['avg_lfc'], "{cond}-lo_{mark}_bs100.log2fc.avg.bw"),
    cond=sorted(CONDITIONS),
    mark=sorted(MARKS)
)
COND_PAIRS = list(combinations(sorted(CONDITIONS),2))
COND_COMB = [f'{a}-vs-{b}' for a,b in COND_PAIRS]
CHROMATIN_STATES = ["actProm", "poiProm", "actEnh", "poiEnh"]

##############################################################################
# RULE ALL
##############################################################################

rule all:
    input:
        expand(
            os.path.join(DIR_DICT['deseq2'],
                         f"{{mark}}_{REPJOIN}_{{cond_comb}}.allPeaks"),
            mark=MARKS,
            cond_comb=COND_COMB
        ),
        expand(
            os.path.join(DIR_DICT['heatmap'],
                         "{cond}_" + REPJOIN + "_{state}.heatmap.pdf"),
            cond=CONDITIONS_LO,
            state=CHROMATIN_STATES
        ),
        expand(
            os.path.join(DIR_DICT['heatmap'],
                         f"{{cond_comb}}_{REPJOIN}_geneBodies.heatmap.pdf"),
            cond_comb=COND_COMB
        ),
        expand(
            os.path.join(DIR_DICT['heatmap'],
                         f"H3K4me3_{{cond}}_{REPJOIN}_geneBodies.hm.gz"),
            cond=CONDITIONS_LO)

rule clean:
    params:
        count_dir=DIR_DICT['counts'],
        merge_dir=DIR_DICT['mergePeaks']
    shell:
        '''
        rm -r {params.count_dir}
        rm -r {params.merge_dir}
        '''

##############################################################################
# Get peaks and counts
##############################################################################

#rule consensus_peaks:
#    input:
#        expand(os.path.join(MACS3_DIR,
#                            "{rep}-{cond}-lo_{mark}_peaks.broadPeak"),
#               rep=REPS,
#               cond=CONDITIONS,
#               mark="{mark}")
#    output:
#        saf=os.path.join(DIR_DICT['mergePeaks'],
#                         f"{{mark}}_{REPJOIN}_consensusPeaks.saf"),
#        cond_beds=expand(
#            os.path.join(DIR_DICT['mergePeaks'],
#                         "{mark}_{cond_lo}_" + REPJOIN + "_mergePeaks.bed"),
#            cond_lo=CONDITIONS_LO,
#            allow_missing=True
#        )
#    params:
#        dir_out=DIR_DICT['mergePeaks'],
#        min_overlap=MIN_OVERLAP
#    log:
#        "logs/consensus_peaks_{mark}.log"
#    shell:
#        '''
#        python3 scripts/01_consensusPeaks.py \
#            {input} \
#            {params.dir_out} \
#            {params.min_overlap} \
#            > {log} 2>&1
#        '''

rule consensus_peaks:
    input:
        expand(os.path.join(MACS3_DIR,
                            "{rep}-{cond}-lo_{mark}_peaks.broadPeak"),
               rep=REPS,
               cond=CONDITIONS,
               mark="{mark}")
    output:
        bed=os.path.join(DIR_DICT['mergePeaks'],
                         f"{{mark}}_{REPJOIN}_consensusPeaks.bed"),
        saf=os.path.join(DIR_DICT['mergePeaks'],
                         f"{{mark}}_{REPJOIN}_consensusPeaks.saf")
    params:
        dir_out=DIR_DICT['mergePeaks']
    log:
        "logs/consensus_peaks_{mark}.log"
    shell:
        '''
        cat {input} \
            | sort -k1,1 -k2,2n \
            | bedtools merge \
            > {output.bed} 2> {log}

        awk 'OFS="\t" {{print $1"_"$2"_"$3, $1, $2, $3, "."}}' {output.bed} \
            > {output.saf} 2>> {log}
        '''

rule feature_count:
    input:
        bam_file=os.path.join(BAM_DIR,
                              f"{{rep}}-{{cond}}-lo_{{mark}}{BAM_SUFFIX}"),
        saf_file=os.path.join(DIR_DICT['mergePeaks'],
                              f"{{mark}}_{REPJOIN}_consensusPeaks.saf")
    output:
        os.path.join(DIR_DICT['counts'], 
                     f"{{rep}}-{{cond}}_{{mark}}{BAM_SUFFIX.replace('.bam','.counts')}")
    log:
        "logs/featurecount_{rep}-{cond}_{mark}.log"
    resources:
        mem_mb=16000,
        runtime=60
    shell:
        '''
        python3 scripts/02_featureCount.py \
            {input.bam_file} \
            {input.saf_file} \
            {output} \
            > {log} 2>&1
        '''

rule count_matrix:
    input:
        count_files=expand(
            os.path.join(DIR_DICT['counts'],
                         "{rep}-{cond}_{mark}"+BAM_SUFFIX.replace('.bam', '.counts')),
            rep=REPS,
            cond=CONDITIONS,
            mark="{mark}",
            allow_missing=True
        )
    output:
        os.path.join(DIR_DICT['countMatrix'], 
                     f"{{mark}}_{REPJOIN}_consensusPeaks.countMatrix")
    log:
        "logs/count_matrix_{mark}.log"
    shell:
        '''
        python3 scripts/03_countMatrix.py \
            {input.count_files} \
            {output} \
            > {log} 2>&1
        '''

##############################################################################
# Run DESeq2
##############################################################################

rule deseq2:
    input:
        os.path.join(DIR_DICT['countMatrix'], 
                     f"{{mark}}_{REPJOIN}_consensusPeaks.countMatrix")
    output:
        os.path.join(DIR_DICT['deseq2'], 
                     f"{{mark}}_{REPJOIN}_{{cond_comb}}.allPeaks")
    wildcard_constraints:
        cond_comb="|".join(COND_COMB)
    params:
        alpha=ALPHA,
        lfc_threshold=LFC_THRESHOLD,
        dir_out=DIR_DICT['deseq2']
    log:
        "logs/deseq2_{mark}_{cond_comb}.log"
    resources:
        mem_mb=64000,
        runtime=480
    shell:
        '''
        python3 scripts/04_deseq2.py \
            {input} \
            {params.alpha} \
            {params.lfc_threshold} \
            {params.dir_out} \
            > {log} 2>&1
        '''

##############################################################################
# Average bigwigs
##############################################################################

rule average_thor:
    input:
        expand(os.path.join(DIR_DICT['thor'],
                "THOR_{{mark}}_{{cond_comb}}-{{samp}}-{reps}.bw"),
                reps=THOR_REPS)
    output:
        os.path.join(DIR_DICT['avg_thor'],
                    "THOR_{mark}_{cond_comb}_{samp}_bs100.avg.bw")
    log:
        'logs/average_thor_{mark}_{cond_comb}_{samp}.log'
    resources:
        mem_mb=128000,
        runtime=240
    shell:
        '''
        bigwigAverage \
            -b {input} \
            -bs 100 \
            -p max \
            -o {output} \
            > {log} 2>&1
        '''

rule average_lfc:
    input:
        expand(os.path.join(LFC_BW,
                "{reps}-{{cond}}-lo_{{mark}}_bs100.log2fc.bigwig"),
                reps=REPS)
    output:
        os.path.join(DIR_DICT['avg_lfc'],
                    "{cond}-lo_{mark}_bs100.log2fc.avg.bw")
    log:
        'logs/average_lfc_{cond}_{mark}.log'
    resources:
        mem_mb=128000,
        runtime=240
    shell:
        '''
        bigwigAverage \
            -b {input} \
            -bs 100 \
            -p max \
            -o {output} \
            > {log} 2>&1
        '''
        
##############################################################################
# Run deeptools comparisons
##############################################################################

rule chromatin_state_actProm:
    input:
        h3k4me3=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K4me3_{{cond}}_{REPJOIN}_mergePeaks.bed"),
        h3k27ac=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K27ac_{{cond}}_{REPJOIN}_mergePeaks.bed")
    output:
        os.path.join(DIR_DICT['intersect'], f"{{cond}}_{REPJOIN}_actProm.bed")
    wildcard_constraints:
        cond="|".join(CONDITIONS_LO)
    log:
        "logs/chromatin_state_actProm_{cond}.log"
    shell:
        '''
        bedtools intersect \
            -a {input.h3k4me3} \
            -b {input.h3k27ac} \
            -wa \
            > {output} \
            2> {log}
        '''

rule intersect_h3k4me3_gene_bodies:
    input:
        peaks=os.path.join(DIR_DICT['intersect'], f"{{cond}}_{REPJOIN}_actProm.bed"),
        gene_bodies=GENE_BODIES_BED
    output:
        os.path.join(DIR_DICT['intersect'],
                     f"H3K4me3_{REPJOIN}_{{cond}}.allPeaks.geneBodies.bed")
    wildcard_constraints:
        cond_comb="|".join(COND_COMB)
    log:
        "logs/intersect_h3k4me3_gene_bodies_{cond}.log"
    shell:
        '''
        bedtools intersect \
            -a {input.gene_bodies} \
            -b {input.peaks} \
            -wa \
            > {output} \
            2> {log}
        '''

rule chromatin_state_poiProm:
    input:
        h3k4me3=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K4me3_{{cond}}_{REPJOIN}_mergePeaks.bed"),
        h3k27ac=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K27ac_{{cond}}_{REPJOIN}_mergePeaks.bed")
    output:
        os.path.join(DIR_DICT['intersect'], f"{{cond}}_{REPJOIN}_poiProm.bed")
    wildcard_constraints:
        cond="|".join(CONDITIONS_LO)
    log:
        "logs/chromatin_state_poiProm_{cond}.log"
    shell:
        '''
        bedtools intersect \
            -a {input.h3k4me3} \
            -b {input.h3k27ac} \
            -v \
            > {output} \
            2> {log}
        '''

rule chromatin_state_actEnh:
    input:
        h3k4me1=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K4me1_{{cond}}_{REPJOIN}_mergePeaks.bed"),
        h3k27ac=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K27ac_{{cond}}_{REPJOIN}_mergePeaks.bed")
    output:
        os.path.join(DIR_DICT['intersect'], f"{{cond}}_{REPJOIN}_actEnh.bed")
    wildcard_constraints:
        cond="|".join(CONDITIONS_LO)
    log:
        "logs/chromatin_state_actEnh_{cond}.log"
    shell:
        '''
        bedtools intersect \
            -a {input.h3k4me1} \
            -b {input.h3k27ac} \
            -wa \
            > {output} \
            2> {log}
        '''

rule chromatin_state_poiEnh:
    input:
        h3k4me1=os.path.join(DIR_DICT['mergePeaks'],
	                     f"H3K4me1_{{cond}}_{REPJOIN}_mergePeaks.bed"),
        h3k27ac=os.path.join(DIR_DICT['mergePeaks'],
                             f"H3K27ac_{{cond}}_{REPJOIN}_mergePeaks.bed")
    output:
        os.path.join(DIR_DICT['intersect'], f"{{cond}}_{REPJOIN}_poiEnh.bed")
    wildcard_constraints:
        cond="|".join(CONDITIONS_LO)
    log:
        "logs/chromatin_state_poiEnh_{cond}.log"
    shell:
        '''
        bedtools intersect \
            -a {input.h3k4me1} \
            -b {input.h3k27ac} \
            -v \
            > {output} \
            2> {log}
        '''

##############################################################################
# deeptools heatmaps
##############################################################################

rule compute_matrix_chromatin_states:
    input:
        bigwigs=BIGWIGS,
        regions=os.path.join(DIR_DICT['intersect'],
                             "{cond}_" + REPJOIN + "_{state}.bed")
    output:
        os.path.join(DIR_DICT['heatmap'],
                     "{cond}_" + REPJOIN + "_{state}.matrix.gz")
    wildcard_constraints:
        cond="|".join(CONDITIONS_LO),
        state="|".join(CHROMATIN_STATES)
    params:
        b=DT_B,
        a=DT_A,
        bs=DT_BS
    log:
        "logs/compute_matrix_{cond}_{state}.log"
    resources:
        mem_mb=32000,
        runtime=120
    shell:
        '''
        computeMatrix reference-point \
            -S {input.bigwigs} \
            -R {input.regions} \
            -b {params.b} \
            -a {params.a} \
            -bs {params.bs} \
	    -p max \
            -o {output} \
            > {log} 2>&1
        '''

rule plot_heatmap_chromatin_states:
    input:
        os.path.join(DIR_DICT['heatmap'],
                     "{cond}_" + REPJOIN + "_{state}.matrix.gz")
    output:
        heatmap=os.path.join(DIR_DICT['heatmap'],
                             "{cond}_" + REPJOIN + "_{state}.heatmap.pdf"),
        tab=os.path.join(DIR_DICT['heatmap'],
                         "{cond}_" + REPJOIN + "_{state}.hm.gz")
    wildcard_constraints:
        cond="|".join(CONDITIONS_LO),
        state="|".join(CHROMATIN_STATES)
    params:
        labels=" ".join([
            os.path.basename(f).split("_bs100")[0].replace("-lo_", "_")
            for f in BIGWIGS
        ])
    log:
        "logs/plot_heatmap_{cond}_{state}.log"
    shell:
        '''
        plotHeatmap \
            -m {input} \
            -o {output.heatmap} \
	    --colorMap coolwarm \
            --outFileNameMatrix {output.tab} \
            --samplesLabel {params.labels} \
            > {log} 2>&1
        '''

rule compute_matrix_gene_bodies:
    input:
        bigwigs=BIGWIGS,
        regions=GENE_BODIES_BED
    output:
        os.path.join(DIR_DICT['heatmap'],
                     f"{{cond_comb}}_{REPJOIN}_geneBodies.matrix.gz")
    wildcard_constraints:
        cond_comb="|".join(COND_COMB)
    params:
        b=DT_B,
        a=DT_A,
        m=DT_M,
        bs=DT_BS
    log:
        "logs/compute_matrix_gene_bodies_{cond_comb}.log"
    resources:
        mem_mb=32000,
        runtime=120
    shell:
        '''
        computeMatrix scale-regions \
            -S {input.bigwigs} \
            -R {input.regions} \
            -b {params.b} \
            -a {params.a} \
            -m {params.m} \
	    -p max \
            -bs {params.bs} \
            -o {output} \
            > {log} 2>&1
        '''

rule plot_heatmap_gene_bodies:
    input:
        os.path.join(DIR_DICT['heatmap'],
                     f"{{cond_comb}}_{REPJOIN}_geneBodies.matrix.gz")
    output:
        heatmap=os.path.join(DIR_DICT['heatmap'],
                             f"{{cond_comb}}_{REPJOIN}_geneBodies.heatmap.pdf"),
        tab=os.path.join(DIR_DICT['heatmap'],
                         f"{{cond_comb}}_{REPJOIN}_geneBodies.hm.gz")
    wildcard_constraints:
        cond_comb="|".join(COND_COMB)
    params:
        labels=" ".join([
            os.path.basename(f).split("_bs100")[0].replace("-lo_", "_")
            for f in BIGWIGS
        ])
    log:
        "logs/plot_heatmap_gene_bodies_{cond_comb}.log"
    shell:
        '''
        plotHeatmap \
            -m {input} \
            -o {output.heatmap} \
	    --colorMap coolwarm \
            --outFileNameMatrix {output.tab} \
            --samplesLabel {params.labels} \
            > {log} 2>&1
        '''

rule compute_matrix_k4me3_gene_bodies:
    input:
        bigwigs=BIGWIGS,
        regions=os.path.join(DIR_DICT['intersect'],
                     f"H3K4me3_{REPJOIN}_{{cond}}.allPeaks.geneBodies.bed")
    output:
        os.path.join(DIR_DICT['heatmap'],
                     f"H3K4me3_{{cond}}_{REPJOIN}_geneBodies.matrix.gz")
    wildcard_constraints:
        cond_comb="|".join(COND_COMB)
    params:
        b=DT_B,
        a=DT_A,
        m=DT_M,
        bs=DT_BS
    log:
        "logs/compute_matrix_h3k4me3_gene_bodies_{cond}.log"
    resources:
        mem_mb=32000,
        runtime=120
    shell:
        '''
        computeMatrix scale-regions \
            -S {input.bigwigs} \
            -R {input.regions} \
            -b {params.b} \
            -a {params.a} \
            -m {params.m} \
            -p max \
            -bs {params.bs} \
            -o {output} \
            > {log} 2>&1
        '''

rule plot_heatmap_k4me3_gene_bodies:
    input:
        os.path.join(DIR_DICT['heatmap'],
                     f"H3K4me3_{{cond}}_{REPJOIN}_geneBodies.matrix.gz")
    output:
        heatmap=os.path.join(DIR_DICT['heatmap'],
                             f"H3K4me3_{{cond}}_{REPJOIN}_geneBodies.heatmap.pdf"),
        tab=os.path.join(DIR_DICT['heatmap'],
                         f"H3K4me3_{{cond}}_{REPJOIN}_geneBodies.hm.gz")
    wildcard_constraints:
        cond_comb="|".join(COND_COMB)
    params:
        labels=" ".join([
            os.path.basename(f).split("_bs100")[0].replace("-lo_", "_")
            for f in BIGWIGS
        ])
    log:
        "logs/plot_heatmap_k4me3_gene_bodies_{cond}.log"
    shell:
        '''
        plotHeatmap \
            -m {input} \
            -o {output.heatmap} \
            --colorMap coolwarm \
            --outFileNameMatrix {output.tab} \
            --samplesLabel {params.labels} \
            > {log} 2>&1
        '''

