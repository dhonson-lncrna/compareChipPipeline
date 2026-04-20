import os
import sys

import pandas as pd
import numpy as np

import pybedtools
from collections import defaultdict

def load_peaks(peak_file):
    """
    Load MACS3 narrowPeak or broadPeak file
    """
    cols = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
            'signalValue', 'pValue', 'qValue', 'peak']
    df = pd.read_csv(peak_file, sep='\t', header=None, 
                     names=cols[:min(10, len(open(peak_file).readline().split('\t')))])
    return pybedtools.BedTool.from_dataframe(df[['chrom', 'start', 'end']]).sort()

def rep_shared_peaks(rep_peaks,
                     outname,
                     min_overlap=0.5):
    '''
    Identify shared peaks between replicates.
    
    Parameters:
    -----------
    rep_peaks : list
        List of paths to peak files.
    min_overlap : float
        Minimum reciprocal overlap fraction (default 0.5 = 50%)
    
    Returns:
    --------
    BedTool : Merged peaks present in both replicates
    '''
    # Load peaks
    rep_peaks = [load_peaks(f) for f in rep_peaks]

    # If two replicates, do simple intersection. Otherwise iterate.
    if len(rep_peaks) == 2:
        intersect = rep_peaks[0].intersect(rep_peaks[1], f=min_overlap, r=True)
    else:
        intersect = rep_peaks[0].intersect(rep_peaks[1], f=min_overlap, r=True)
        for i in rep_peaks[2:]:
            intersect = intersect.intersect(i, f=min_overlap, r=True)

    # Merge overlapping intervals
    merged = intersect.sort().merge()
    if outname:
        merged.saveas(outname)

    return merged

def write_saf(bedfile):
    '''
    Writes a simplified annotation format (SAF) file for use with 
    featureCounts
    '''
    saf_out = bedfile.replace('.bed','.saf')
    peaks = pd.read_csv(bedfile,
                        sep='\t',
                        header=None,
                        names=['chrom','start','end'])

    with open(saf_out,'w') as saf:
        saf.write('GeneID\tChr\tStart\tEnd\tStrand\n')
        for i, row in peaks.iterrows():
            saf.write(f"peak_{i}\t{row['chrom']}\t{row['start']+1}\t{row['end']}\t.\n")

    return None

def consensus_peakset(peak_files,
                      dir_out,
                      min_overlap=0.5):
    """
    Create consensus peak set across all samples. In Snakemake,
    implemented as taking only peakfiles for a single target.
    
    Parameters:
    -----------
    dir_in : str
        Directory containing MACS3 narrowPeaks
    dir_out: str
        Directory to output merge files
    conditions : list
        List of condition names
    min_overlap : float
        Minimum overlap for calling shared peaks
    
    Returns:
    --------
    BedTool : Consensus peak set (union of all peaks merged)
    """
    # Get files and experiment info
    narrow_peaks = [os.path.basename(i) for i in peak_files]
    print(narrow_peaks)
    reps = [i.split('-')[0] for i in narrow_peaks]
    repjoin = '-'.join(sorted(set(reps)))
    conds = [i.split('-',1)[1].split('_')[0] for i in narrow_peaks]
    targets = [i.split('_')[1] for i in narrow_peaks]
    
    expt_ls = list(zip(targets,conds,peak_files))
    expt_dict = defaultdict(lambda: defaultdict(list))
    for m,c,f in expt_ls:
        expt_dict[m][c].append(f)

    # Get shared peaks
    all_peaks = {t:[] for t in targets}
    
    for ka, va in expt_dict.items():
        for kb, vb in va.items():
            outname = f'{ka}_{kb}_{repjoin}_mergePeaks.bed'
            outname = os.path.join(dir_out,outname)
            rep_shared_peaks(vb,
                             outname,
                             min_overlap)
            all_peaks[ka].append(outname)
    
    # Combine all condition-specific peaks
    for k, v in all_peaks.items():
        # Combine peaks
        bedtools_ls = [pybedtools.BedTool(fpath) for fpath in v]
        combined = bedtools_ls[0].cat(*bedtools_ls[1:], postmerge=False)
        
        # Sort and merge
        consensus = combined.sort().merge()
        bed_out = f'{k}_{repjoin}_consensusPeaks.bed'
        bed_out = os.path.join(dir_out, bed_out)
        consensus.saveas(bed_out)

        # Write SAF file for featureCounts
        write_saf(bed_out)
    
    return consensus

if __name__ == '__main__':
    args = sys.argv
    if len(args) >= 4:
        consensus_peakset(args[1:-2],
                          args[-2],
                          args[-1])
    else:
        print("Usage: 01_consensusPeaks.py <peak_files> <dir_out> <min_overlap>", file=sys.stderr)
        sys.exit(1)
    
