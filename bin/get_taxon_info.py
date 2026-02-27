import os
import argparse
from Bio import Entrez

cache_dir = os.path.join(os.getenv('TMPDIR', '/tmp'), 'biopython_cache')
os.makedirs(cache_dir, exist_ok=True)
Entrez.local_cache = cache_dir

def get_taxon_info(taxon_id):
    """Fetch taxonomic information from NCBI Taxonomy database."""
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error fetching data: {e}")
        return None

    if not records:
        return None

    record = records[0]
    taxon_info = {
        'TxID': taxon_id,
        'Species': record.get('ScientificName', 'N/A'),
        'Phylum': 'N/A',
        'Class': 'N/A',
        'Order': 'N/A'
    }

    lineage = record.get('LineageEx', [])

    for entry in lineage:
        if entry['Rank'] == 'phylum':
            taxon_info['Phylum'] = entry.get('ScientificName', 'N/A')
        elif entry['Rank'] == 'class':
            taxon_info['Class'] = entry.get('ScientificName', 'N/A')
        elif entry['Rank'] == 'order':
            taxon_info['Order'] = entry.get('ScientificName', 'N/A')

    return taxon_info

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Retrieve taxonomic information for a given NCBI taxon ID')
    parser.add_argument('--taxon_id', help='NCBI Taxonomy ID (e.g., 9606 for Homo sapiens)')
    parser.add_argument('--email', help="Email for querying NCBI database")
    parser.add_argument('--output', help='Output file name', default="Taxon_info.txt")
    args = parser.parse_args()

    Entrez.email = args.email
    result = get_taxon_info(args.taxon_id)
    
    with open(args.output, 'w') as f:
        if result:
            f.write(f"TxID\t{result['TxID']}\n")
            f.write(f"Species\t{result['Species']}\n")
            f.write(f"Phylum\t{result['Phylum']}\n")
            f.write(f"Class\t{result['Class']}\n")
            f.write(f"Order\t{result['Order']}\n")
        else:
            f.write(f"No taxonomic information found for ID: {args.taxon_id}\n")