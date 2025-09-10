import sys
import pickle
from Bio.Seq import Seq
from Bio import SeqIO

def translate_dna(dna_seq, genetic_code):
    seq = Seq(dna_seq.upper())
    return str(seq.translate(table=genetic_code))

def analyze_introns(input_file, genetic_code):
    stats = {
        'total_introns': 0,
        'short_intron': 0,
        'long_intron': 0,
        'short_intron_counts': {0: 0, 1: 0, 2: 0},
        'long_intron_counts': {0: 0, 1: 0, 2: 0},
        'short_intron_counts_with_stop': {0: 0, 1: 0, 2: 0},
        'short_intron_counts_without_stop': {0: 0, 1: 0, 2: 0},
        'long_intron_counts_with_stop': {0: 0, 1: 0, 2: 0},
        'long_intron_counts_without_stop': {0: 0, 1: 0, 2: 0},
        'plotting_data': {
            'with_stop': {0: [], 1: [], 2: []},
            'without_stop': {0: [], 1: [], 2: []}
        }
    }
    
    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq)
        seq_len = len(seq_str)
        mod3 = seq_len % 3
        stats['total_introns'] += 1
        
        protein_seq = translate_dna(seq_str, genetic_code)
        has_stop = '*' in protein_seq
        
        stop_key = 'with_stop' if has_stop else 'without_stop'
        stats['plotting_data'][stop_key][mod3].append(seq_len)
        
        if seq_len <= 120:
            stats['short_intron_counts'][mod3] += 1
            stats['short_intron'] += 1
            stats_key = f'short_intron_counts_{stop_key}'
        else:
            stats['long_intron_counts'][mod3] += 1
            stats['long_intron'] += 1
            stats_key = f'long_intron_counts_{stop_key}'
        stats[stats_key][mod3] += 1
    
    return stats

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_intron_stats.py intron_file genetic_code")
        sys.exit(1)

    intron_file = sys.argv[1]
    genetic_code = int(sys.argv[2])
    
    stats = analyze_introns(intron_file, genetic_code)
    total = stats['total_introns']
    stopless_short_intron = sum(stats["short_intron_counts_without_stop"].values())
    stop_containing_short_intron = sum(stats["short_intron_counts_with_stop"].values())
    
    with open("stop_codon_statistics.tsv", 'w') as f:
        f.write(f"short_intron_3n/short_intron\t{stats["short_intron_counts"][0]} ({stats["short_intron_counts"][0]/stats["short_intron"]*100:.2f}%)\n")
        f.write(f"long_intron_3n/long_intron\t{stats["long_intron_counts"][0]} ({stats["long_intron_counts"][0]/stats["long_intron"]*100:.2f}%)\n")

        f.write(f"stopless_short_intron_3n/stopless_short_intron\t{stats["short_intron_counts_without_stop"][0]} ({stats["short_intron_counts_without_stop"][0]/stopless_short_intron*100:.2f}%)\n")
        f.write(f"stopless_short_intron_3n1/stopless_short_intron\t{stats["short_intron_counts_without_stop"][1]} ({stats["short_intron_counts_without_stop"][1]/stopless_short_intron*100:.2f}%)\n")
        f.write(f"stopless_short_intron_3n2/stopless_short_intron\t{stats["short_intron_counts_without_stop"][2]} ({stats["short_intron_counts_without_stop"][2]/stopless_short_intron*100:.2f}%)\n")
  
        f.write(f"stop_containing_short_intron_3n/short_containing_intron\t{stats['short_intron_counts_with_stop'][0]} ({stats["short_intron_counts_with_stop"][0]/stop_containing_short_intron*100:.2f}%)\n")
        f.write(f"stop_containing_short_intron_3n1/short_containing_intron\t{stats['short_intron_counts_with_stop'][1]} ({stats["short_intron_counts_with_stop"][1]/stop_containing_short_intron*100:.2f}%)\n")
        f.write(f"stop_containing_short_intron_3n2/short_containing_intron\t{stats['short_intron_counts_with_stop'][2]} ({stats["short_intron_counts_with_stop"][2]/stop_containing_short_intron*100:.2f}%)\n")

    # Save plotting data
    with open("intron_plot_data.pkl", 'wb') as pkl_file:
        pickle.dump(stats['plotting_data'], pkl_file)

if __name__ == '__main__':
    main()