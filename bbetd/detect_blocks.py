import csv
from collections import defaultdict

def detect_tandem_blocks(order_tsv: str, diamond_tsv: str, out_blocks_tsv: str, max_gap: int = 5):
    """
    Clusters BLAST hits into connected tandem arrays/blocks on the same chromosome.
    """
    gene_order = {}
    with open(order_tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene_order[row['gene_id']] = (row['chrom'], int(row['order']))

    # Graph adjacency for collinear/nearby hits
    adj = defaultdict(set)
    with open(diamond_tsv, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 2:
                continue
            g1, g2 = cols[0].strip(), cols[1].strip()
            if g1 == g2:
                continue
            if g1 in gene_order and g2 in gene_order:
                chr1, ord1 = gene_order[g1]
                chr2, ord2 = gene_order[g2]
                if chr1 == chr2 and abs(ord1 - ord2) - 1 <= max_gap:
                    adj[g1].add(g2)
                    adj[g2].add(g1)

    # Connected components
    visited = set()
    blocks = []
    block_id = 1
    
    for gene in gene_order:
        if gene in visited or gene not in adj:
            continue
        
        # BFS / DFS
        component = []
        queue = [gene]
        visited.add(gene)
        
        while queue:
            curr = queue.pop(0)
            component.append(curr)
            for neighbor in adj[curr]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
                    
        if len(component) > 1:
            # Sort component by genomic order
            component.sort(key=lambda g: gene_order[g][1])
            chrom = gene_order[component[0]][0]
            blocks.append((f"BLOCK_{block_id:05d}", chrom, ",".join(component)))
            block_id += 1

    with open(out_blocks_tsv, 'w') as f:
        f.write("block_id\tchrom\tgenes\n")
        for b_id, ch, glist in blocks:
            f.write(f"{b_id}\t{ch}\t{glist}\n")