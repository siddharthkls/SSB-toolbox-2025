from collections import defaultdict
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import pandas as pd

# Load your 1000 Pfam domains
background_df = pd.read_csv("D:/Toolbox/Result.csv")
background_pfam_ids = set(background_df['Pfam-ThreeThousand Control'].dropna().astype(str))

def load_pfam2go(filename):
    pfam2go = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('Pfam:'):
                try:
                    parts = line.strip().split('>')
                    pfam_part = parts[0].split()[0]       # 'Pfam:PF00004'
                    pfam_acc = pfam_part.split(':')[1]     # 'PF00004'
                    
                    go_desc = parts[1].split(';')[0].strip()  # 'GO:ATP binding'
                    go_id = parts[1].split(';')[1].strip()    # 'GO:0005524'
                    full_go = f"{go_id} {go_desc}"

                    if pfam_acc not in pfam2go:
                        pfam2go[pfam_acc] = []
                    pfam2go[pfam_acc].append(full_go)
                except Exception as e:
                    print(f"[Warning] Skipped line due to error: {e}")
    return pfam2go


def read_pfam_sectors(file_path):
    sectors = []
    current = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                if current:
                    sectors.append(current)
                    current = []
            else:
                current.append(line)
        if current:
            sectors.append(current)
    return sectors


def map_sectors_to_go(sectors, pfam2go_map):
    results = []
    for i, sector in enumerate(sectors, 1):
        go_terms = []
        for pfam_id in sector:
            go_terms.extend(pfam2go_map.get(pfam_id, []))
        results.append((f"sector {i}", sorted(set(go_terms))))
    return results


def compute_enrichment(go_sectors, pfam_sectors, pfam2go_map, go2pfam_map, background_pfam_ids):
    N = len(background_pfam_ids)  # Restrict background to 1000 Pfams
    all_results = []

    for i, (sector_name, go_terms) in enumerate(go_sectors):
        sector_pfams = set(pfam_sectors[i]) & background_pfam_ids  # Use only Pfams from 1000
        n = len(sector_pfams)
        term_stats = []

        for go_term in go_terms:
            pfams_with_go = go2pfam_map.get(go_term, set()) & background_pfam_ids
            K = len(pfams_with_go)
            k = len(sector_pfams & pfams_with_go)

            # Hypergeometric p-value
            pval = hypergeom.sf(k - 1, N, K, n)
            term_stats.append((go_term, k, K, pval))

        pvals = [x[3] for x in term_stats]
        qvals = multipletests(pvals, method='fdr_bh')[1] if pvals else []

        enriched = []
        for j, (go_term, k, K, pval) in enumerate(term_stats):
            qval = qvals[j]
            enriched.append((go_term, k, K, pval, qval))
        all_results.append((sector_name, enriched))
    
    return all_results



# === MAIN SCRIPT ===

# Load mappings and sector data
pfam2go_map = load_pfam2go("D:/Toolbox/Pfam2go.txt")
pfam_sectors = read_pfam_sectors("D:/Toolbox/Pfam - ThreeThousand.txt")
go_sectors = map_sectors_to_go(pfam_sectors, pfam2go_map)

# Reverse map: GO term → list of Pfam domains
go2pfam_map = defaultdict(set)
for pfam, go_list in pfam2go_map.items():
    for go in go_list:
        go2pfam_map[go].add(pfam)

# Enrichment analysis
enrichment_results = compute_enrichment(go_sectors, pfam_sectors, pfam2go_map, go2pfam_map, background_pfam_ids)

# Output
output_file = "D:/Toolbox/Pfam_GO_3000.txt"
with open(output_file, 'w') as out:
    for sector_name, enriched_terms in enrichment_results:
        out.write(f"{sector_name}:\n")
        if enriched_terms:
            for go, k, K, pval, qval in sorted(enriched_terms, key=lambda x: x[4]):  # sort by q-value
                status = "ENRICHED" if qval < 0.05 else ""
                out.write(f"  - {go} | in sector: {k}, in background: {K}, p={pval:.4g}, q={qval:.4g} {status}\n")
        else:
            out.write("  No GO terms found\n")
        out.write("\n")

print(f"[Done] Enrichment analysis written to: {output_file}")

import csv
import re

# Define input and output filenames
input_filename = 'D:/Toolbox/Pfam_GO_3000.txt'  # Replace with your actual filename
output_filename = 'Pfam_GO_3000.csv'

# Regular expressions to extract sector name and GO term description with q-value
sector_pattern = re.compile(r"^(sector \d+):")
go_pattern = re.compile(r"(GO:\d{7}) ([^|]+)\s+\|.*q=([0-9.]+) .*ENRICHED")

extracted_data = []
current_sector = None

with open(input_filename, 'r') as infile:
    for line in infile:
        line = line.strip()
        # Check if line is a sector name
        sector_match = sector_pattern.match(line)
        if sector_match:
            current_sector = sector_match.group(1)
            continue
        
        # Check if line contains enriched GO term info
        go_match = go_pattern.search(line)
        if go_match and current_sector:
           go_id = go_match.group(1)
           description = go_match.group(2).strip()
           q_value = float(go_match.group(3))
           extracted_data.append((current_sector, go_id, description, q_value))


# Write to CSV with sector column
with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['sector', 'GO_ID', 'GO Term Description', 'FDR (q-value)'])
    writer.writerows(extracted_data)

print(f"Extracted {len(extracted_data)} enriched GO terms with sectors into '{output_filename}'")


import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from goatools.obo_parser import GODag

# Load updated CSV with GO_ID column
df = pd.read_csv("Pfam_GO_3000.csv")

# Load the OBO DAG from the downloaded file
from goatools.obo_parser import GODag
go_dag = GODag("D:/Toolbox/go-basic.obo")  # use local path

# Map each GO_ID to its namespace
def get_namespace(go_id):
    try:
        term = go_dag[go_id]
        if term.namespace == 'biological_process':
            return 'BP'
        elif term.namespace == 'molecular_function':
            return 'MF'
        elif term.namespace == 'cellular_component':
            return 'CC'
        else:
            return 'Other'
    except KeyError:
        return 'Unknown'

df['Namespace'] = df['GO_ID'].apply(get_namespace)

# Save annotated file
df.to_csv("Pfam_GO_3000_annotated.csv", index=False)
print("[✓] GO terms classified and saved to 'Pfam_GO_3000_annotated.csv'")
