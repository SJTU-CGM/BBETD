import shutil
import subprocess

def run_diamond_blastp(pep_fasta: str, db_prefix: str, hits_tsv: str, 
                       threads: int = 8, evalue: float = 1e-2, 
                       matrix: str = "BLOSUM45", diamond_bin: str = "diamond"):
    """Build Diamond DB and run self-blastp search."""
    if not shutil.which(diamond_bin):
        raise FileNotFoundError(f"Diamond binary '{diamond_bin}' not found in PATH.")

    # 1. Build database
    makedb_cmd = [diamond_bin, "makedb", "--in", pep_fasta, "-d", db_prefix]
    subprocess.run(makedb_cmd, check=True)

    # 2. Run blastp
    blast_cmd = [
        diamond_bin, "blastp",
        "-d", f"{db_prefix}.dmnd",
        "-q", pep_fasta,
        "--ultra-sensitive",
        "--matrix", matrix,
        "--evalue", str(evalue),
        "--comp-based-stats", "0",
        "--max-target-seqs", "0",
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qlen", "slen",
        "-o", hits_tsv,
        "--threads", str(threads)
    ]
    subprocess.run(blast_cmd, check=True)