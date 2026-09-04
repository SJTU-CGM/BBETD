import os
import re

def clean_id(gid: str) -> str:
    gid = gid.strip()
    gid = re.sub(r'^(gene:|transcript:|rna-|maker-|cds:|exon:)', '', gid)
    return gid.strip()

def classify_tandem_pairs(order_tsv: str, diamond_tsv: str, blocks_tsv: str, 
                          out_pairs: str, max_gap: int = 5, 
                          max_evalue: float = 1e-2, min_max_cov: float = 20.0):
    """
    Applies Dual-Rule Classification:
      - Rule 1: Immediately adjacent (gap == 0) -> Accept if BLAST hit passes coverage/e-value.
      - Rule 2: Intervening gap (1 <= gap <= max_gap) -> Must belong to a detected tandem block.
    """
    # 1. Load Gene Orders
    gene_order = {}
    with open(order_tsv, 'r') as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                gid = clean_id(parts[0])
                chrom = parts[1].strip()
                ord_idx = int(parts[4])
                gene_order[gid] = (chrom, ord_idx)

    # 2. Filter Valid Diamond Hits
    valid_hits = set()
    with open(diamond_tsv, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 8:
                a, b = clean_id(cols[0]), clean_id(cols[1])
                if a == b:
                    continue
                
                length = float(cols[3])
                evalue = float(cols[4])
                qlen = float(cols[6])
                slen = float(cols[7])
                
                if qlen <= 0 or slen <= 0:
                    continue
                
                qcov = (length / qlen) * 100.0
                scov = (length / slen) * 100.0
                max_cov = max(qcov, scov)
                
                if evalue <= max_evalue and max_cov >= min_max_cov:
                    valid_hits.add((min(a, b), max(a, b)))

    # 3. Load Tandem Blocks
    block_pairs = set()
    if os.path.exists(blocks_tsv):
        with open(blocks_tsv, 'r') as f:
            f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    genes = [clean_id(g) for g in parts[2].split(',') if g.strip()]
                    for i in range(len(genes)):
                        for j in range(i + 1, len(genes)):
                            if genes[i] != genes[j]:
                                block_pairs.add((min(genes[i], genes[j]), max(genes[i], genes[j])))

    # 4. Classify Pairs
    predicted_pairs = set()
    for (g1, g2) in valid_hits:
        if g1 in gene_order and g2 in gene_order:
            chrom1, ord1 = gene_order[g1]
            chrom2, ord2 = gene_order[g2]
            if chrom1 == chrom2:
                gap = abs(ord1 - ord2) - 1
                pair_key = (min(g1, g2), max(g1, g2))
                
                # Rule 1: Adjacent
                if gap == 0:
                    predicted_pairs.add((g1, g2))
                # Rule 2: In-block with gap <= max_gap
                elif 1 <= gap <= max_gap:
                    if pair_key in block_pairs:
                        predicted_pairs.add((g1, g2))

    with open(out_pairs, 'w') as f:
        f.write('Duplicate 1\tDuplicate 2\n')
        for g1, g2 in sorted(predicted_pairs):
            f.write(f"{g1}\t{g2}\n")

    return len(predicted_pairs)