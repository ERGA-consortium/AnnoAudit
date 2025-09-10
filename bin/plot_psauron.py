import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_fasta(file_path):
    """Parse FASTA file and return dictionary of description: protein_length"""
    desc_to_length = {}
    current_desc = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_desc is not None:
                    seq = ''.join(current_seq).replace('*', '')
                    desc_to_length[current_desc] = len(seq)
                current_desc = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_desc and current_seq:
            seq = ''.join(current_seq).replace('*', '')
            desc_to_length[current_desc] = len(seq)
    return desc_to_length

def main():
    parser = argparse.ArgumentParser(description='Generate protein score distribution with length analysis')
    parser.add_argument('--input', required=True, help='Path to the input CSV file')
    parser.add_argument('--fasta', required=True, help='Path to the protein FASTA file')
    parser.add_argument('--output', '-o', default='PSAURON_plot', help='Base name for output files')
    args = parser.parse_args()

    df = pd.read_csv(args.input, skiprows=2)
    desc_length = parse_fasta(args.fasta)
    df['protein_length'] = df['description'].map(desc_length)

    bins = np.arange(0, 1.05, 0.05)
    df['score_bin'] = pd.cut(df['in-frame_score'], bins=bins, include_lowest=True)
    df['bin_mid'] = df['score_bin'].apply(lambda x: x.mid)

    count_data = df.groupby(['bin_mid', 'psauron_is_protein'], observed=False).size().unstack(fill_value=0)
    length_data = df.groupby('bin_mid', observed=False)['protein_length'].mean().reset_index()

    sns.set_theme(style="whitegrid")
    fig, ax1 = plt.subplots(figsize=(12, 7))

    ax1.bar(bins[:-1]+0.025, count_data[True], width=0.045, color='#1f77b4', label='True')
    ax1.bar(bins[:-1]+0.025, count_data[False], width=0.045, color='#ff7f0e', 
            bottom=count_data[True], label='False')

    ax1.set_xlabel('In-Frame Score Bins', fontsize=12)
    ax1.set_ylabel('Number of Proteins', fontsize=12)
    ax1.set_xticks(bins[:-1] + 0.025)
    ax1.set_xticklabels([f"{b:.2f}-{b+0.05:.2f}" for b in bins[:-1]], rotation=45, ha='right')
    ax1.set_xlim(0, 1.0)

    ax2 = ax1.twinx()
    ax2.plot(length_data['bin_mid'], length_data['protein_length'], color='purple',
            marker='D', linestyle='-', label='Mean Protein Length')

    ax2.set_ylabel('Mean Protein Length (aa)', fontsize=12)
    ax2.grid(False)

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, 
              bbox_to_anchor=(1.15, 1), loc='upper left')

    plt.title('Protein Distribution and Length Relationship by In-Frame Score', 
             pad=20, fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.savefig(f'{args.output}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{args.output}.pdf', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()