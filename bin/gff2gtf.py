import sys
import pandas as pd

def parse_attributes(attribute_str):
    """Parse GFF attribute string into a dictionary"""
    attr_dict = {}
    for attribute in attribute_str.split(';'):
        attribute = attribute.strip()
        if attribute and '=' in attribute:
            key, value = attribute.split('=', 1)
            attr_dict[key.strip()] = value.strip()
    return attr_dict

def parse_gff(gff_file):
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv(gff_file, sep='\t', header=None, names=columns, comment='#')
    return df

def convert_gff_to_gtf(gff_df):
    gtf_lines = []
    transcript_to_gene = {}
    
    # First pass: build transcript to gene mapping
    for _, row in gff_df.iterrows():
        attr_dict = parse_attributes(row['attribute'])
        feature = row['feature']
        
        if feature in ['mRNA', 'transcript']:
            transcript_id = attr_dict.get('ID', '')
            gene_id = attr_dict.get('Parent', '')
            
            # Sometimes Parent might not be set, check gene attribute
            if not gene_id:
                gene_id = attr_dict.get('gene', '')
                if gene_id:
                    # Format gene_id properly
                    gene_id = f'gene-{gene_id}' if not gene_id.startswith('gene-') else gene_id
            
            if transcript_id and gene_id:
                transcript_to_gene[transcript_id] = gene_id
    
    # Second pass: generate GTF lines
    for _, row in gff_df.iterrows():
        attr_dict = parse_attributes(row['attribute'])
        feature = row['feature']
        
        if feature == 'gene':
            continue
            
        elif feature in ['mRNA', 'transcript']:
            transcript_id = attr_dict.get('ID', '')
            gene_id = transcript_to_gene.get(transcript_id, '')
            
            if not gene_id:
                gene_id = attr_dict.get('gene', '')
                if gene_id:
                    gene_id = f'gene-{gene_id}' if not gene_id.startswith('gene-') else gene_id
            
            if not (transcript_id and gene_id):
                continue
            
            attributes = [
                f'gene_id "{gene_id}"',
                f'transcript_id "{transcript_id}"'
            ]
            attr_str = '; '.join(attributes) + ';'
            gtf_line = (
                f"{row.seqname}\t{row.source}\ttranscript\t"
                f"{row.start}\t{row.end}\t{row.score}\t"
                f"{row.strand}\t{row.frame}\t{attr_str}"
            )
            gtf_lines.append(gtf_line)
            
        elif feature in ['exon', 'CDS']:
            parent_transcript = attr_dict.get('Parent', '')
            gene_id = transcript_to_gene.get(parent_transcript, '')
            
            # If we can't find gene_id from transcript, try other methods
            if not gene_id:
                # Try from gene attribute
                gene_id = attr_dict.get('gene', '')
                if gene_id:
                    gene_id = f'gene-{gene_id}' if not gene_id.startswith('gene-') else gene_id
                # Try from transcript_to_gene using different transcript ID formats
                elif parent_transcript:
                    # Check if parent_transcript is in dictionary
                    for key in transcript_to_gene:
                        if key in parent_transcript or parent_transcript in key:
                            gene_id = transcript_to_gene[key]
                            break
            
            if not gene_id:
                # Skip if we can't find gene_id
                continue
            
            attributes = [f'gene_id "{gene_id}"']
            
            if feature == 'exon':
                exon_id = attr_dict.get('exon_id', attr_dict.get('ID', ''))
                if exon_id:
                    attributes.append(f'exon_id "{exon_id}"')
            elif feature == 'CDS':
                cds_id = attr_dict.get('ID', '')
                if cds_id:
                    attributes.append(f'cds_id "{cds_id}"')
            
            if parent_transcript:
                attributes.append(f'transcript_id "{parent_transcript}"')
            
            attr_str = '; '.join(attributes) + ';'
            gtf_line = (
                f"{row.seqname}\t{row.source}\t{feature.lower()}\t"
                f"{row.start}\t{row.end}\t{row.score}\t"
                f"{row.strand}\t{row.frame}\t{attr_str}"
            )
            gtf_lines.append(gtf_line)
    
    return gtf_lines

def write_gtf_file(gtf_lines, output_file):
    with open(output_file, 'w') as f:
        for line in gtf_lines:
            f.write(line + '\n')

def main(gff_file, output_file):
    gff_df = parse_gff(gff_file)
    gtf_lines = convert_gff_to_gtf(gff_df)
    write_gtf_file(gtf_lines, output_file)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python gff_to_gtf.py input.gff output.gtf")
        sys.exit(1)

    gff_file = sys.argv[1]
    output_file = sys.argv[2]
    main(gff_file, output_file)