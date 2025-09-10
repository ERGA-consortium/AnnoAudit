import gzip
import argparse
import pandas as pd
from Bio import SeqIO

def parse_gff(gff_file):
    """Parse GFF3 file and return DataFrames for genes, mRNAs, exons, CDS, introns, and transcript lengths"""
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    
    opener = gzip.open if gff_file.endswith('.gz') else open
    with opener(gff_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t', comment='#', header=None, names=columns)
    
    # Convert frame column: replace '.' with '0' and convert to integer
    df['frame'] = df['frame'].replace('.', '0')
    df['frame'] = pd.to_numeric(df['frame'], errors='coerce').fillna(0).astype(int)
    
    df['ID'] = df['attribute'].str.extract(r'ID=([^;]+)')
    df['Parent'] = df['attribute'].str.extract(r'Parent=([^;]+)')
    
    mrnas = df[df['feature'] == 'mRNA'].copy()
    transcript_to_gene = dict(zip(mrnas['ID'], mrnas['Parent']))
    
    df['gene_id'] = df['Parent'].map(transcript_to_gene)
    df.loc[df['feature'] == 'gene', 'gene_id'] = df.loc[df['feature'] == 'gene', 'ID']
    
    genes = df[df['feature'] == 'gene'].copy()
    exons = df[df['feature'] == 'exon'].copy()
    cds = df[df['feature'] == 'CDS'].copy()
    
    # Assign transcript_id to exons and CDS
    exons['transcript_id'] = exons['Parent']
    cds['transcript_id'] = cds['Parent']
    
    # Calculate introns and transcript lengths
    exons_sorted = exons.sort_values(['transcript_id', 'start'])
    introns = []
    transcript_lengths = {}
    
    for transcript_id, group in exons_sorted.groupby('transcript_id'):
        exons_list = group.to_dict('records')
        merged_exons = merge_exons(exons_list)
        if len(merged_exons) < 1:
            continue
        
        # Calculate transcript length (sum of merged exons)
        transcript_length = sum(e['end'] - e['start'] + 1 for e in merged_exons)
        transcript_lengths[transcript_id] = transcript_length
        
        # Get CDS features for this transcript
        cds_for_transcript = cds[cds['transcript_id'] == transcript_id].to_dict('records')
        strand = exons_list[0]['strand'] if exons_list else '+'
        
        # For each merged exon, find CDS features within it and compute donor phase
        merged_exons_with_phase = []
        for exon in merged_exons:
            cds_in_exon = []
            for cds_feat in cds_for_transcript:
                if (cds_feat['start'] >= exon['start'] and 
                    cds_feat['end'] <= exon['end']):
                    cds_in_exon.append(cds_feat)
            
            # Sort CDS features by genomic order (strand-aware)
            if strand == '+':
                cds_in_exon.sort(key=lambda x: x['start'])
            else:
                cds_in_exon.sort(key=lambda x: x['start'], reverse=True)
            
            # Compute cumulative phase at end of exon
            current_phase = None
            for feat in cds_in_exon:
                length = feat['end'] - feat['start'] + 1
                if current_phase is None:
                    current_phase = feat['frame']
                current_phase = (current_phase + length) % 3
            
            # If no CDS, mark phase as '.' (non-coding)
            donor_phase = current_phase if current_phase is not None else '.'
            merged_exons_with_phase.append({
                'start': exon['start'],
                'end': exon['end'],
                'donor_phase': donor_phase
            })
        
        # Create introns between consecutive merged exons
        chrom = merged_exons[0]['seqname'] if merged_exons else None
        gene_id = merged_exons[0]['gene_id'] if merged_exons else None
        for i in range(1, len(merged_exons_with_phase)):
            prev = merged_exons_with_phase[i-1]
            curr = merged_exons_with_phase[i]
            intron_start = prev['end'] + 1
            intron_end = curr['start'] - 1
            if intron_start > intron_end:
                continue
            
            introns.append({
                'seqname': chrom,
                'feature': 'intron',
                'start': intron_start,
                'end': intron_end,
                'strand': strand,
                'gene_id': gene_id,
                'transcript_id': transcript_id,
                'attribute': f'Parent={transcript_id}',
                'donor_phase': prev['donor_phase']
            })
    
    introns_df = pd.DataFrame(introns)
    return genes, mrnas, exons, cds, introns_df, transcript_lengths

def merge_exons(exons):
    """Merge overlapping exons for a transcript."""
    if not exons:
        return []
    sorted_exons = sorted(exons, key=lambda x: x['start'])
    merged = [sorted_exons[0].copy()]
    for exon in sorted_exons[1:]:
        last = merged[-1]
        if exon['start'] <= last['end']:
            merged[-1]['end'] = max(last['end'], exon['end'])
        else:
            merged.append(exon.copy())
    return merged

def merge_overlapping_introns(introns):
    """Merge overlapping or adjacent introns."""
    if not introns:
        return []
    sorted_introns = sorted(introns, key=lambda x: (x[0], x[1], x[2]))
    merged = [list(sorted_introns[0])]
    for intron in sorted_introns[1:]:
        last = merged[-1]
        if intron[0] == last[0] and intron[1] <= last[2]:
            merged[-1][2] = max(last[2], intron[2])
        else:
            merged.append(list(intron))
    return merged

def calculate_statistics(genes, mrnas, exons, cds, introns, genome_size):
    stats = {}
    
    # Gene-level statistics
    stats['num_genes'] = len(genes)
    gene_lengths = genes['end'] - genes['start'] + 1
    stats['mean_gene_length'] = round(gene_lengths.mean(), 2)
    stats['median_gene_length'] = round(gene_lengths.median(), 2)
    
    # Transcript-level statistics
    stats['num_transcripts'] = len(mrnas)
    transcripts_per_gene = mrnas.groupby('Parent').size()
    stats['mean_transcripts_per_gene'] = round(transcripts_per_gene.mean(), 2)
    stats['median_transcripts_per_gene'] = round(transcripts_per_gene.median(), 2)
    
    # Exon statistics
    total_exons = len(exons)
    all_transcript_ids = mrnas['ID']
    exon_counts = exons.groupby('transcript_id').size().reindex(all_transcript_ids, fill_value=0)
    multi_exon_counts = exon_counts[exon_counts >= 2]
    stats['num_exons'] = total_exons
    stats['mean_exons_per_transcript'] = round(exon_counts.mean(), 2)
    stats['median_exons_per_transcript'] = round(exon_counts.median(), 2)
    stats['mean_exons_per_multiexon_transcript'] = round(multi_exon_counts.mean(), 2) if not multi_exon_counts.empty else 0.0
    stats['median_exons_per_multiexon_transcript'] = round(multi_exon_counts.median(), 2) if not multi_exon_counts.empty else 0.0
    num_transcripts_no_introns = (exon_counts == 1).sum()
    stats['num_transcripts_without_introns'] = f"{num_transcripts_no_introns} ({num_transcripts_no_introns/len(exon_counts)*100:.2f}%)"
    
    # Exon mod3 statistics
    exons['exon_length'] = exons['end'] - exons['start'] + 1
    exons['mod3'] = exons['exon_length'] % 3
    mod3_counts = exons['mod3'].value_counts()
    stats['num_exon_3n'] = f"{mod3_counts.get(0, 0)} ({mod3_counts.get(0, 0)/total_exons*100:.2f}%)"
    stats['num_exon_3n1'] = f"{mod3_counts.get(1, 0)} ({mod3_counts.get(1, 0)/total_exons*100:.2f}%)"
    stats['num_exon_3n2'] = f"{mod3_counts.get(2, 0)} ({mod3_counts.get(2, 0)/total_exons*100:.2f}%)"
    
    # CDS statistics
    cds['cds_length'] = cds['end'] - cds['start'] + 1
    cds_per_transcript = cds.groupby('transcript_id')['cds_length'].sum().reindex(all_transcript_ids, fill_value=0)
    stats['mean_cds_length'] = round(cds_per_transcript.mean(), 2)
    stats['median_cds_length'] = round(cds_per_transcript.median(), 2)
    total_cds = cds_per_transcript.sum()
    stats['total_cds_length'] = total_cds
    stats['percentage_cds_coverage'] = f"{total_cds/genome_size*100:.2f}%"
    
    # Intron statistics
    intron_counts = introns.groupby('transcript_id').size().reindex(all_transcript_ids, fill_value=0)
    stats['num_introns'] = len(introns)
    intron_lengths = introns['end'] - introns['start'] + 1
    stats['mean_intron_length'] = round(intron_lengths.mean(), 2)
    stats['median_intron_length'] = round(intron_lengths.median(), 2)
    stats['mean_introns_per_transcript'] = round(intron_counts.mean(), 2)
    stats['median_introns_per_transcript'] = round(intron_counts.median(), 2)
    
    return stats

def extract_intron_sequences(introns_df, fasta_file, output_file):
    with open(fasta_file, "r") as handle:
        genome = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    
    with open(output_file, 'w') as out_f:
        for _, row in introns_df.iterrows():
            chrom = row['seqname']
            start = row['start']
            end = row['end']
            strand = row['strand']
            gene_id = row['gene_id']
            transcript_id = row['transcript_id']
            donor_phase = row['donor_phase']
            
            try:
                seq = genome[chrom].seq[start-1:end]
                if strand == '-':
                    seq = seq.reverse_complement()
                # Include donor_phase in the header
                out_f.write(f">{gene_id}|{transcript_id}|{chrom}:{start}-{end}({strand})|{donor_phase}\n{seq}\n")
            except KeyError:
                continue

def generate_bed_files(introns_df, mrnas, transcript_lengths):
    # Prepare all introns for isoforms_introns.bed
    all_introns_bed = [
        (row['seqname'], row['start'] - 1, row['end'])
        for _, row in introns_df.iterrows()
    ]
    merged_all = merge_overlapping_introns(all_introns_bed)
    with open('isoforms_introns.bed', 'w') as f:
        for chrom, start, end in merged_all:
            f.write(f"{chrom}\t{start}\t{end}\n")
    
    # Prepare longest isoform introns
    transcript_info = pd.DataFrame({
        'transcript_id': list(transcript_lengths.keys()),
        'length': list(transcript_lengths.values())
    })
    mrnas.rename(columns={'ID': 'transcript_id'}, inplace=True)
    transcript_info = transcript_info.merge(
        mrnas[['transcript_id', 'Parent']], on='transcript_id', how='left'
    )
    transcript_info.rename(columns={'Parent': 'gene_id'}, inplace=True)
    
    longest_transcripts = transcript_info.loc[
        transcript_info.groupby('gene_id')['length'].idxmax()
    ]
    longest_transcript_ids = longest_transcripts['transcript_id'].tolist()
    
    longest_introns = introns_df[
        introns_df['transcript_id'].isin(longest_transcript_ids)
    ]
    longest_introns_bed = [
        (row['seqname'], row['start'] - 1, row['end'])
        for _, row in longest_introns.iterrows()
    ]
    merged_longest = merge_overlapping_introns(longest_introns_bed)
    with open('longest_isoform_introns.bed', 'w') as f:
        for chrom, start, end in merged_longest:
            f.write(f"{chrom}\t{start}\t{end}\n")

def main(gff_file, fasta_file):
    genes, mrnas, exons, cds, introns, transcript_lengths = parse_gff(gff_file)
    
    with open(fasta_file, "r") as handle:
        genome_size = sum(len(rec.seq) for rec in SeqIO.parse(handle, "fasta"))
    
    stats = calculate_statistics(genes, mrnas, exons, cds, introns, genome_size)
    pd.Series(stats).to_csv("statistics.tsv", sep='\t', header=False)
    
    extract_intron_sequences(introns, fasta_file, "intron_sequences.fasta")
    generate_bed_files(introns, mrnas, transcript_lengths)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GFF analysis, intron extraction, and BED generation')
    parser.add_argument('gff_file', help='Input GFF3 file')
    parser.add_argument('fasta_file', help='Input genome FASTA file')
    args = parser.parse_args()
    
    main(args.gff_file, args.fasta_file)