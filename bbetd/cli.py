import argparse
import sys
import os
import time
from bbetd.prep import run_gffread, process_pep_and_order
from bbetd.align import run_diamond_blastp
from bbetd.detect_blocks import detect_tandem_blocks
from bbetd.classify import classify_tandem_pairs
from bbetd import __version__

def main():
    parser = argparse.ArgumentParser(
        description="BBETD: Block-Based Extended Tandem Duplication Gene Identifier",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-g", "--genome", required=True, help="Path to genome FASTA file")
    parser.add_argument("-a", "--gff", required=True, help="Path to annotation GFF3/GTF file")
    parser.add_argument("-o", "--outdir", default="./bbetd_results", help="Output directory")
    parser.add_argument("-p", "--prefix", default="bbetd", help="Prefix for output filenames")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of CPU threads for Diamond")
    
    # Diamond & Filtering Options
    parser.add_argument("--evalue", type=float, default=1e-2, help="E-value threshold for Diamond blastp")
    parser.add_argument("--min-cov", type=float, default=20.0, help="Minimum query/subject max-coverage percentage")
    parser.add_argument("--max-gap", type=int, default=5, help="Maximum number of intervening genes (k <= MAX_GAP)")
    parser.add_argument("--matrix", default="BLOSUM45", help="Diamond alignment scoring matrix")
    
    # Executable Paths
    parser.add_argument("--diamond-bin", default="diamond", help="Path to diamond executable")
    parser.add_argument("--gffread-bin", default="gffread", help="Path to gffread executable")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    start_time = time.time()

    raw_pep = os.path.join(args.outdir, f"{args.prefix}_raw.pep")
    clean_pep = os.path.join(args.outdir, f"{args.prefix}.pTpG.pep")
    order_tsv = os.path.join(args.outdir, f"{args.prefix}.Order.tsv")
    dmnd_db = os.path.join(args.outdir, f"{args.prefix}.pTpG_db")
    hits_tsv = os.path.join(args.outdir, f"{args.prefix}_diamond_hits.tsv")
    blocks_tsv = os.path.join(args.outdir, f"{args.prefix}_tandem_blocks.tsv")
    final_pairs = os.path.join(args.outdir, f"{args.prefix}_bbetd_tdg.pairs")

    print("[1/4] Extracting proteins & calculating gene orders...")
    run_gffread(args.gff, args.genome, raw_pep, args.gffread_bin)
    process_pep_and_order(args.gff, raw_pep, clean_pep, order_tsv)

    print("[2/4] Running Diamond self-blastp alignment...")
    run_diamond_blastp(
        pep_fasta=clean_pep,
        db_prefix=dmnd_db,
        hits_tsv=hits_tsv,
        threads=args.threads,
        evalue=args.evalue,
        matrix=args.matrix,
        diamond_bin=args.diamond_bin
    )

    print("[3/4] Detecting tandem duplication blocks...")
    detect_tandem_blocks(order_tsv, hits_tsv, blocks_tsv, max_gap=args.max_gap)

    print("[4/4] Performing Dual-Rule classification...")
    num_pairs = classify_tandem_pairs(
        order_tsv=order_tsv,
        diamond_tsv=hits_tsv,
        blocks_tsv=blocks_tsv,
        out_pairs=final_pairs,
        max_gap=args.max_gap,
        max_evalue=args.evalue,
        min_max_cov=args.min_cov
    )

    elapsed = time.time() - start_time
    print("==========================================================================")
    print(f"[✓] Completed in {elapsed:.2f}s!")
    print(f"[✓] High-confidence tandem duplicate pairs found: {num_pairs}")
    print(f"[✓] Final results saved to: {final_pairs}")
    print("==========================================================================")

if __name__ == "__main__":
    main()