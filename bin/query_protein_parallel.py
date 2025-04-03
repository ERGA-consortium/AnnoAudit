import argparse
import asyncio
import aiohttp
from Bio import Entrez, SeqIO
from io import StringIO, BytesIO
import json
import hashlib
import math

def calculate_bloom_filter_size(n: int, p: float, k: int) -> int:
    if p <= 0 or p >= 1:
        raise ValueError("Probability p must be between 0 and 1 (exclusive).")
    m = - (n * k) / math.log(1 - (p ** (1 / k)))
    return math.ceil(m)

class BloomFilter:
    def __init__(self, size: int, hash_count: int):
        self.size = size
        self.bit_array = [0] * size
        self.hash_count = hash_count

    def _hashes(self, item: str):
        for i in range(self.hash_count):
            hash_result = hashlib.md5(f"{item}{i}".encode()).hexdigest()
            yield int(hash_result, 16) % self.size

    def add(self, item: str):
        for hash_value in self._hashes(item):
            self.bit_array[hash_value] = 1

    def __contains__(self, item: str) -> bool:
        return all(self.bit_array[hash_value] for hash_value in self._hashes(item))

def extract_species(description: str) -> str:
    """Extract species name from FASTA header"""
    if '[' in description and ']' in description:
        parts = description.rsplit('[', 1)
        if len(parts) > 1:
            species_part = parts[-1].split(']', 1)[0]
            return species_part.strip()
    return "Unknown species"

async def fetch_with_retries(session, url, retries=3, delay=1):
    for i in range(retries):
        try:
            async with session.get(url) as response:
                response.raise_for_status()
                content = await response.text()
                await response.release()
                return content
        except aiohttp.ClientError as e:
            if i < retries - 1:
                await asyncio.sleep(delay)
            else:
                print(f"Failed to fetch {url}: {e}")
                return None

async def fetch_protein_sequences(email, taxon_id, target_count=100000, batch_size=6000, 
                                chunk_size=300, max_concurrent_requests=10, 
                                bloom_filter_error_rate=0.005, num_hash=5):
    Entrez.email = email
    sequences = []
    current_taxon_id = str(taxon_id)
    species_counts = {}

    n = target_count
    p = bloom_filter_error_rate
    k = num_hash
    bloom_filter_size = calculate_bloom_filter_size(n, p, k)
    fetched_ids = BloomFilter(size=bloom_filter_size, hash_count=k)

    connector = aiohttp.TCPConnector(
        limit_per_host=max_concurrent_requests,
        force_close=True
    )
    
    async with aiohttp.ClientSession(connector=connector) as session:
        while len(sequences) < target_count and current_taxon_id is not None:
            # Fetch taxonomy information for display
            taxonomy_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={current_taxon_id}&retmode=xml"
            taxonomy_results = await fetch_with_retries(session, taxonomy_url)
            
            if not taxonomy_results:
                print(f"Failed to fetch taxonomy for {current_taxon_id}")
                break

            try:
                taxonomy_data = Entrez.read(BytesIO(taxonomy_results.encode('utf-8')))
                taxon_name = taxonomy_data[0]['ScientificName']
                lineage = taxonomy_data[0]['LineageEx']
            except Exception as e:
                print(f"Error parsing taxonomy data for {current_taxon_id}: {e}")
                break

            print(f"\nProcessing {taxon_name} (TaxID: {current_taxon_id})")

            total_fetched = 0
            while True:
                search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                            f"db=protein&term=txid{current_taxon_id}[Organism:exp]"
                            f"&retstart={total_fetched}&retmax={batch_size}&retmode=json")
                
                search_results = await fetch_with_retries(session, search_url)
                if not search_results:
                    break

                try:
                    search_data = json.loads(search_results)
                    protein_ids = search_data.get('esearchresult', {}).get('idlist', [])
                except Exception as e:
                    print(f"Error parsing search results: {e}")
                    break

                if not protein_ids:
                    break

                new_ids = [pid for pid in protein_ids if pid not in fetched_ids]
                if not new_ids:
                    break

                id_chunks = [new_ids[i:i + chunk_size] for i in range(0, len(new_ids), chunk_size)]
                chunk_tasks = []

                for chunk in id_chunks:
                    if len(chunk_tasks) >= max_concurrent_requests:
                        results = await asyncio.gather(*chunk_tasks)
                        for fetch_results in results:
                            if fetch_results:
                                sequences_batch = list(SeqIO.parse(StringIO(fetch_results), "fasta"))
                                for seq_record in sequences_batch:
                                    if seq_record.id not in fetched_ids:
                                        sequences.append(seq_record)
                                        fetched_ids.add(seq_record.id)
                                        species_name = extract_species(seq_record.description)
                                        species_counts[species_name] = species_counts.get(species_name, 0) + 1
                        chunk_tasks = []

                    fetch_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                                 f"db=protein&id={','.join(chunk)}&rettype=fasta&retmode=text")
                    chunk_tasks.append(fetch_with_retries(session, fetch_url))

                if chunk_tasks:
                    results = await asyncio.gather(*chunk_tasks)
                    for fetch_results in results:
                        if fetch_results:
                            sequences_batch = list(SeqIO.parse(StringIO(fetch_results), "fasta"))
                            for seq_record in sequences_batch:
                                if seq_record.id not in fetched_ids:
                                    sequences.append(seq_record)
                                    fetched_ids.add(seq_record.id)
                                    species_name = extract_species(seq_record.description)
                                    species_counts[species_name] = species_counts.get(species_name, 0) + 1

                total_fetched += len(protein_ids)
                if len(sequences) >= target_count:
                    break

            # Move up to parent taxon
            parent_id = None
            if lineage and len(lineage) > 1:
                parent_id = lineage[-2]['TaxId']
                print(f"Moving up to parent taxon: {lineage[-2]['ScientificName']} ({parent_id})")
            current_taxon_id = parent_id

    await connector.close()
    return sequences[:target_count], species_counts

def main():
    parser = argparse.ArgumentParser(description='Fetch protein sequences from NCBI with species counts')
    parser.add_argument("-e", '--email', type=str, required=True, help='Email address for NCBI Entrez')
    parser.add_argument("-t", '--taxon_id', type=int, required=True, help='Starting taxon ID')
    parser.add_argument("-c", '--target_count', type=int, default=100000, help='Target sequence count')
    parser.add_argument("-b", '--batch_size', type=int, default=6000, help='Search batch size')
    parser.add_argument("-k", '--chunk_size', type=int, default=300, help='Fetch chunk size')
    parser.add_argument("-m", '--max_concurrent_requests', type=int, default=10, help='Max concurrent requests')
    parser.add_argument("-n", '--num_hash', type=int, default=5, help='Bloom filter hash functions')
    parser.add_argument("-r", '--error_rate', type=float, default=0.005, help='Bloom filter error rate')
    
    args = parser.parse_args()

    sequences, species_report = asyncio.run(
        fetch_protein_sequences(
            args.email, args.taxon_id, 
            args.target_count, args.batch_size,
            args.chunk_size, args.max_concurrent_requests,
            args.error_rate, args.num_hash
        )
    )

    with open('NCBI_protein_sequences.fasta', 'w') as f:
        SeqIO.write(sequences, f, 'fasta')

    with open('species_counts.tsv', 'w') as f:
        f.write("Species\tCount\n")
        for species, count in species_report.items():
            f.write(f"{species}\t{count}\n")

    print(f"\nTotal sequences retrieved: {len(sequences)}")
    print("Species counts saved to species_counts.tsv")

if __name__ == "__main__":
    main()