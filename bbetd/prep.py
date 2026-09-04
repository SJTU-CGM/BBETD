import os
import re
import shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

def clean_gene_id(id_str: str) -> str:
    """Clean gene/transcript IDs from common GFF prefixes and isoforms."""
    id_m = re.search(r'ID=([^;]+)', id_str)
    if id_m:
        gid = id_m.group(1).strip()
    else:
        name_m = re.search(r'Name=([^;]+)', id_str)
        gid = name_m.group(1).strip() if name_m else id_str.strip()
    
    gid = re.sub(r'^(gene:|transcript:|rna-|maker-|cds:|exon:)', '', gid)
    gid = re.sub(r'\.\d{1,2}_dup$', '_dup', gid)
    gid = re.sub(r'\.\d{1,2}$', '', gid)
    return gid.strip()

def run_gffread(gff_path: str, genome_fa: str, raw_pep_out: str, gffread_bin: str = "gffread"):
    """Extract protein sequences using gffread."""
    if not shutil.which(gffread_bin):
        raise FileNotFoundError(f"gffread binary '{gffread_bin}' not found in PATH.")

    cmd = [gffread_bin, "-W", gff_path, "-g", genome_fa, "-y", raw_pep_out]
    subprocess.run(cmd, check=True)

def process_pep_and_order(gff_path: str, raw_pep_path: str, out_pep: str, out_order: str):
    """
    1. Filter the longest isoform per gene.
    2. Build the chromosome-based Gene Order table.
    """
    # A. Select longest protein sequence per gene
    gene_to_longest = {}
    for rec in SeqIO.parse(raw_pep_path, 'fasta'):
        clean_seq = str(rec.seq).replace('*', '').replace('.', '')
        if not clean_seq:
            continue
            
        gid = clean_gene_id(rec.id)
        seq_len = len(clean_seq)
        
        if gid not in gene_to_longest or seq_len > gene_to_longest[gid][0]:
            rec.id = gid
            rec.name = ''
            rec.description = ''
            rec.seq = Seq(clean_seq)
            gene_to_longest[gid] = (seq_len, rec)

    longest_records = [val[1] for val in gene_to_longest.values()]
    SeqIO.write(longest_records, out_pep, 'fasta')

    # B. Generate Gene Order
    order_lines = ['gene_id\tchrom\tstart\tend\torder\n']
    chrom_order = {}

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2].lower() != 'gene':
                continue
            chrom, start, end, attr = parts[0], int(parts[3]), int(parts[4]), parts[8]
            gid = clean_gene_id(attr)
            
            chrom_order[chrom] = chrom_order.get(chrom, 0) + 1
            order_lines.append(f"{gid}\t{chrom}\t{start}\t{end}\t{chrom_order[chrom]}\n")

    with open(out_order, 'w') as f:
        f.writelines(order_lines)