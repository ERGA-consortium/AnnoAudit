import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt

def plot_distribution(data, title, ax, max_length=500):
    mod_classes = [0, 1, 2]
    colors = ['#0f0f0f', '#ff7f0e', '#2ca02c']
    labels = ['3n', '3n+1', '3n+2']
    
    all_lengths = []
    for mc in mod_classes:
        all_lengths.extend(data[mc])
    
    if not all_lengths:
        ax.set_title(f"{title}\n(No data)")
        return
    
    bins = np.arange(0, min(max(all_lengths), max_length) + 1, 10)
    
    for mc, color, label in zip(mod_classes, colors, labels):
        lengths = data[mc]
        if not lengths:
            continue
        counts, _ = np.histogram(lengths, bins=bins)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        ax.plot(bin_centers, counts, label=label, 
                color=color, linewidth=2, marker='o', markersize=4)
    
    ax.set_title(title, fontsize=12)
    ax.set_xlabel('Intron Length (nt)', fontsize=10)
    ax.set_ylabel('Count', fontsize=10)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_length)

def main():
    parser = argparse.ArgumentParser(description='Plot intron length distributions')
    parser.add_argument('--input_pkl', help='Pickle file from extract_intron_stats.py')
    parser.add_argument('--max_length', type=int, default=150,
                        help='Maximum length to display (default: 150)')
    args = parser.parse_args()

    # Load precomputed data
    with open(args.input_pkl, 'rb') as pkl_file:
        plot_data = pickle.load(pkl_file)
    
    # Create plots
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    plot_distribution(
        plot_data['without_stop'], 
        "Introns Without Stop Codons", 
        axs[0], 
        args.max_length
    )
    plot_distribution(
        plot_data['with_stop'], 
        "Introns Containing Stop Codons", 
        axs[1], 
        args.max_length
    )
    
    plt.tight_layout()
    plt.savefig('intron_distributions.pdf')
    plt.savefig('intron_distributions.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    main()