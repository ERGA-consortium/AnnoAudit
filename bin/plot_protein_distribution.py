import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse

def main(input_file):
    # Read and process data
    qlen, slen = [], []
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            qlen.append(int(parts[2]))
            slen.append(int(parts[3]))
    
    df = pd.DataFrame({
        'Length': qlen + slen,
        'Proteome': ['Reference']*len(qlen) + ['Predicted']*len(slen)
    })

    # Set global style
    sns.set_theme(style="whitegrid", font_scale=1.2)
    plt.rcParams['figure.dpi'] = 300
    palette = {'Reference': '#1f77b4', 'Predicted': '#ff7f0e'}

    # Linear-scale plot
    plt.figure(figsize=(12, 8))
    ax = sns.histplot(
        data=df,
        x='Length',
        hue='Proteome',
        kde=True,
        bins=100,
        common_norm=False,
        palette=palette,
        alpha=0.3,
        linewidth=0.5,
        legend=True
    )
    ax.legend_.set_title('Proteome')
    ax.legend_.set_frame_on(True)
    plt.title('Protein Length Distribution (Linear Scale)', pad=20)
    plt.xlabel('Protein Length (amino acids)')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig('Protein_distribution.pdf')
    plt.savefig('Protein_distribution.png')
    plt.close()

    # Log-scale plot
    plt.figure(figsize=(12, 8))
    ax = sns.histplot(
        data=df,
        x='Length',
        hue='Proteome',
        kde=True,
        bins=100,
        common_norm=False,
        palette=palette,
        alpha=0.3,
        linewidth=0.5,
        legend=True
    )
    plt.xscale('log')
    ax.legend_.set_title('Proteome')
    ax.legend_.set_frame_on(True)
    plt.title('Protein Length Distribution (Log Scale)', pad=20)
    plt.xlabel('Protein Length (log scale)')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig('Protein_distribution_log.pdf')
    plt.savefig('Protein_distribution_log.png')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot protein length distributions from DIAMOND blastp output.')
    parser.add_argument("-i", "--input_file", required=True, help="Input DIAMOND blastp results file")
    args = parser.parse_args()
    main(args.input_file)