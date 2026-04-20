import os
import argparse
import subprocess


def make_config(BAMFILES,
                INPUTS,
                DESTDIR,
                MARK,
                REP1,
                REP2,
                CHROMSIZES='/resnick/groups/carnegie_poc/mcfallngai-lab/dhonson/genomes/escolopes_v2/fastaRename/escolopes_v2simple.chrom.sizes',
                GENOME='/resnick/groups/carnegie_poc/mcfallngai-lab/dhonson/genomes/escolopes_v2/fastaRename/escolopes_v2simple.fasta'):
    REP1BAM = [f + '\n' for f in BAMFILES if REP1 in os.path.basename(f)]
    REP1INPUT = [f + '\n' for f in INPUTS if REP1 in os.path.basename(f)]
    REP2BAM = [f + '\n' for f in BAMFILES if REP2 in os.path.basename(f)]
    REP2INPUT = [f + '\n' for f in INPUTS if REP2 in os.path.basename(f)]

    FNAME = f'THOR_{MARK}_{REP1}-vs-{REP2}.config'
    FPATH = os.path.join(DESTDIR, FNAME)

    with open(FPATH, 'w') as f:
        f.writelines('#rep1\n')
        f.writelines(REP1BAM)
        f.writelines('#rep2\n')
        f.writelines(REP2BAM)
        f.writelines('#chrom_sizes\n')
        f.writelines(CHROMSIZES + '\n')
        f.writelines('#genome\n')
        f.writelines(GENOME + '\n')
        f.writelines('#inputs1\n')
        f.writelines(REP1INPUT)
        f.writelines('#inputs2\n')
        f.writelines(REP2INPUT)

    return FPATH


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=(
            'Generate a rgt-THOR config file for differential ChIP-seq analysis '
            'and run THOR with standard parameters (--merge, --no-correction, '
            '-p 0.01, --binsize 500).'
        )
    )
    parser.add_argument(
        '--bams',
        nargs='+',
        required=True,
        metavar='BAM',
        dest='BAMFILES',
        help=(
            'Space-separated list of paths to ChIP BAM files (pre-filtered for the '
            'correct histone mark). REP1 and REP2 identifiers are used to assign '
            'files to each condition.'
        )
    )
    parser.add_argument(
        '--inputs',
        nargs='+',
        required=True,
        metavar='INPUT',
        dest='INPUTS',
        help=(
            'Space-separated list of paths to input/control BAM files. '
            'REP1 and REP2 identifiers are used to assign files to each condition.'
        )
    )
    parser.add_argument(
        'DESTDIR',
        help='Output directory where the .config file (and THOR results) will be written.'
    )
    parser.add_argument(
        'MARK',
        help='Histone mark identifier used for naming the output config file (e.g. H3K4me1, H3K27ac).'
    )
    parser.add_argument(
        'REP1',
        help='String identifying replicate 1 in BAM filenames (e.g. a condition name or sample ID).'
    )
    parser.add_argument(
        'REP2',
        help='String identifying replicate 2 in BAM filenames (e.g. a condition name or sample ID).'
    )
    parser.add_argument(
        '--chromsizes',
        default='/resnick/groups/carnegie_poc/mcfallngai-lab/dhonson/genomes/escolopes_v2/fastaRename/escolopes_v2simple.chrom.sizes',
        dest='CHROMSIZES',
        help='Path to chromosome sizes file. Defaults to the E. scolopes v2 chrom.sizes.'
    )
    parser.add_argument(
        '--genome',
        default='/resnick/groups/carnegie_poc/mcfallngai-lab/dhonson/genomes/escolopes_v2/fastaRename/escolopes_v2simple.fasta',
        dest='GENOME',
        help='Path to genome FASTA file. Defaults to the E. scolopes v2 assembly.'
    )

    args = parser.parse_args()

    config_path = make_config(
        args.BAMFILES,
        args.INPUTS,
        args.DESTDIR,
        args.MARK,
        args.REP1,
        args.REP2,
        args.CHROMSIZES,
        args.GENOME,
    )

    # Derive PREFIX by stripping the .config extension from the config filepath
    PREFIX = os.path.splitext(config_path)[0]

    cmd = [
        'rgt-THOR',
        '--name', PREFIX,
        '--merge',
        '--no-correction',
        '-p', '0.01',
        '--binsize', '500',
        '--parallel',
        f'{PREFIX}.config',
    ]
    print(f'Running: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
