import os
import sys
import subprocess

def count_peak_frags(bam_file,
                     saf_file,
                     count_out):
    '''
    Uses featureCounts to extract read counts in peak regions
    specified by a SAF file.

    Parameters:
    -----------
    bam_file : str
        Path to sorted and indexed BAM file.
    saf_file : str
        Path to consensus peakset SAF file.
    dir_out : str
        Directory where output will be deposited.

    Returns:
    --------
    A tab-separated counts file.
    '''
    # Count file
    cmd = ['featureCounts',
           '-p',           # Paired end mode
           '-F', 'SAF',    # SAF input
           '-s', '0',      # Unstranded mode
           '-a', saf_file,
           '-o', count_out,
           bam_file]
    
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            capture_output=True, 
            text=True
        )
    except subprocess.CalledProcessError as e:
        print(f"featureCounts failed with error:\n{e.stderr}")
        raise

    return None

if __name__ == '__main__':
    args = sys.argv
    if len(args) == 4:
        count_peak_frags(args[1],
                         args[2],
                         args[3])
    else:
        print("Usage: 02_featureCount.py <bam_file> <saf_file> <count_out>", file=sys.stderr)
        sys.exit(1)
    