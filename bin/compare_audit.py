import json
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.gridspec import GridSpec
import re

def parse_json_files(json_files):
    """Parse multiple JSON files and extract key metrics for comparison"""
    species_data = []
    
    for file_path in json_files:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        # Parse BUSCO strings
        busco_annot = data['BUSCO']['BUSCO Annot']
        busco_assem = data['BUSCO']['BUSCO Assem']
        
        def parse_busco_string(busco_str):
            """Parse BUSCO string format: C:94.5[S:94.0,D:0.6],F:0.8,M:4.7"""
            try:
                c_match = re.search(r'C:([\d.]+)', busco_str)
                completeness = float(c_match.group(1)) if c_match else np.nan
                
                s_match = re.search(r'S:([\d.]+)', busco_str)
                d_match = re.search(r'D:([\d.]+)', busco_str)
                single_copy = float(s_match.group(1)) if s_match else np.nan
                duplicated = float(d_match.group(1)) if d_match else np.nan
                
                f_match = re.search(r'F:([\d.]+)', busco_str)
                m_match = re.search(r'M:([\d.]+)', busco_str)
                fragmented = float(f_match.group(1)) if f_match else np.nan
                missing = float(m_match.group(1)) if m_match else np.nan
                
                return completeness, single_copy, duplicated, fragmented, missing
            except Exception as e:
                print(f"Error parsing BUSCO string '{busco_str}': {e}")
                return np.nan, np.nan, np.nan, np.nan, np.nan
        
        annot_complete, annot_single, annot_duplicated, annot_fragmented, annot_missing = parse_busco_string(busco_annot)
        assem_complete, assem_single, assem_duplicated, assem_fragmented, assem_missing = parse_busco_string(busco_assem)
        
        def parse_percentage(percent_str):
            """Extract percentage value from string like '12825 (90.35%)'"""
            try:
                return float(percent_str.split('(')[1].split('%')[0])
            except:
                return np.nan
        
        def parse_rnaseq_percentage(percent_str):
            """Extract percentage value from RNA-Seq strings"""
            try:
                if '(' in percent_str:
                    return float(percent_str.split('(')[1].split('%')[0])
                else:
                    return float(percent_str.strip('%'))
            except:
                return np.nan
        
        # Parse intron phase percentages from General Statistics
        def parse_intron_percentage(percent_str):
            """Extract percentage value from intron statistics strings"""
            try:
                return float(percent_str.split(' ')[1].strip('()%'))
            except:
                return np.nan
        
        # Parse Best Reciprocal Hits percentages
        def parse_brh_percentage(percent_str):
            """Extract percentage value from Best Reciprocal Hits strings"""
            try:
                return float(percent_str.split(' ')[1].strip('()%'))
            except:
                return np.nan
        
        species_info = {
            'ToLID': data['Taxon Info']['ToLID'],
            'Species': data['Taxon Info']['Species'],
            'Phylum': data['Taxon Info']['Phylum'],
            'Class': data['Taxon Info']['Class'],
            'Order': data['Taxon Info']['Order'],
            
            # Basic statistics
            'num_genes': int(data['Summary']['num_genes']),
            'num_transcripts': int(data['Summary']['num_transcripts']),
            'transcripts_without_introns_pct': float(data['Summary']['num_transcripts_without_introns'].split(' ')[1].strip('()%')),
            'transcripts_without_introns_count': int(data['Summary']['num_transcripts_without_introns'].split(' ')[0]),
            
            # Structural features - using median values
            'median_exons_per_transcript': float(data['General Statistics']['median_exons_per_transcript']),
            'median_exons_per_multiexon_transcript': float(data['General Statistics']['median_exons_per_multiexon_transcript']),
            'median_cds_length': float(data['Summary']['median_cds_length']),
            'median_intron_length': float(data['Summary']['median_intron_length']),
            
            # BUSCO metrics
            'busco_annotation_completeness': annot_complete,
            'busco_annotation_single': annot_single,
            'busco_annotation_duplicated': annot_duplicated,
            'busco_assembly_completeness': assem_complete,
            'busco_assembly_single': assem_single,
            'busco_assembly_duplicated': assem_duplicated,
            'delta_busco': float(data['BUSCO']['Delta BUSCO'].strip('%')),
            'delta_busco_abs': abs(float(data['BUSCO']['Delta BUSCO'].strip('%'))),  # Absolute value for radar plot
            
            # Intron phase statistics for Plot B
            'stopless_3n': parse_intron_percentage(data['General Statistics']['stopless_short_intron_3n/stopless_short_intron']),
            'stopless_3n1': parse_intron_percentage(data['General Statistics']['stopless_short_intron_3n1/stopless_short_intron']),
            'stopless_3n2': parse_intron_percentage(data['General Statistics']['stopless_short_intron_3n2/stopless_short_intron']),
            'stop_containing_3n': parse_intron_percentage(data['General Statistics']['stop_containing_short_intron_3n/short_containing_intron']),
            'stop_containing_3n1': parse_intron_percentage(data['General Statistics']['stop_containing_short_intron_3n1/short_containing_intron']),
            'stop_containing_3n2': parse_intron_percentage(data['General Statistics']['stop_containing_short_intron_3n2/short_containing_intron']),
            
            # Quality metrics
            'psauron_score': float(data['PSAURON']['psauron_score']),
            
            # OMArk metrics
            'omark_consistent': parse_percentage(data['OMArk']['total_consistent']),
            'omark_missing': parse_percentage(data['OMArk']['missing']),
            
            # RNA-Seq support
            'gene_support': parse_percentage(data['RNASeq']['num_gene_supported']),
            'exon_support': parse_percentage(data['RNASeq']['num_exon_supported']),
            'intron_support': parse_percentage(data['RNASeq']['num_exact_intron_boundary']),
            
            # Best Reciprocal Hits
            'split_genes_pct': parse_brh_percentage(data['Best Reciprocal Hits']['num_split_gene_length_<0.75x']),
            'fusion_genes_pct': parse_brh_percentage(data['Best Reciprocal Hits']['num_fusion_gene_length_>1.25x'])
        }
        species_data.append(species_info)
    
    return pd.DataFrame(species_data)

def create_comprehensive_plot(df):
    """Create one comprehensive figure with all 4 subplots"""
    
    # Create consistent ordering for all plots - sort by ToLID
    df_sorted = df.sort_values('ToLID')
    
    # Create the main figure
    fig = plt.figure(figsize=(22, 18))
    
    # Create grid specification
    gs = GridSpec(2, 2, figure=fig)  # 2 rows, 2 columns
    
    # Plot 1: Radar Plot (Top-left)
    ax1 = fig.add_subplot(gs[0, 0], projection='polar')
    
    # Calculate reverse values for split and fusion genes (higher is better for radar)
    df_sorted['split_genes_pct_reverse'] = 100 - df_sorted['split_genes_pct']
    df_sorted['fusion_genes_pct_reverse'] = 100 - df_sorted['fusion_genes_pct']
    
    # Updated radar plot metrics
    metrics = [
        'omark_consistent',
        'psauron_score', 
        'busco_annotation_completeness',
        'gene_support',
        'exon_support', 
        'intron_support',
        'delta_busco_abs',
        'split_genes_pct_reverse',
        'fusion_genes_pct_reverse'
    ]
    
    metric_labels = [
        'OMArk\nConsistent', 
        'PSAURON\nScore',
        'BUSCO\nComplete', 
        'RNA-Seq Gene\nSupport',
        'RNA-Seq Exon\nSupport',
        'RNA-Seq Intron\nSupport',
        'ΔBUSCO\n(100-|Δ|)',  # Changed to show it's inverted
        'Split Genes\n(100-% split)',
        'Fusion Genes\n(100-% fusion)'
    ]
    
    # Normalize data for radar plot (0-100 scale)
    normalized_data = []
    for i, metric in enumerate(metrics):
        # All metrics are already percentages or scores (0-100), except delta_busco_abs
        if metric == 'delta_busco_abs':
            # Invert delta BUSCO so that smaller delta (better) is higher on radar
            # Convert to 0-100 scale where 0 delta = 100, 100 delta = 0
            normalized = 100 - df_sorted[metric]
        elif metric in ['split_genes_pct_reverse', 'fusion_genes_pct_reverse']:
            # These are already 100 - percentage, so they're already in correct direction
            normalized = df_sorted[metric].clip(0, 100)
        else:
            # Other metrics: clip to 0-100 range
            normalized = df_sorted[metric].clip(0, 100)
        normalized_data.append(normalized)
    
    # Set up radar plot angles
    angles = np.linspace(0, 2*np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]  # Close the circle
    
    # Use a colormap for different species in radar plot
    colors_radar = plt.cm.tab10(np.linspace(0, 1, len(df_sorted)))
    
    for idx, (tolid, color) in enumerate(zip(df_sorted['ToLID'], colors_radar)):
        values = [normalized_data[i].iloc[idx] for i in range(len(metrics))]
        values += values[:1]  # Close the circle
        
        ax1.plot(angles, values, 'o-', linewidth=2, label=tolid, 
                color=color, markersize=4)
        ax1.fill(angles, values, alpha=0.1, color=color)
    
    # Add metric labels
    ax1.set_xticks(angles[:-1])
    ax1.set_xticklabels(metric_labels, fontsize=9)  # Slightly smaller font for more metrics
    ax1.set_ylim(0, 100)
    
    # Add concentric circles
    ax1.set_yticks([20, 40, 60, 80, 100])
    ax1.set_yticklabels(['20', '40', '60', '80', '100'], fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    ax1.set_title('A: Quality Metrics Radar Plot', size=12, fontweight='bold', pad=20)
    
    # Move radar plot legend to below the plot
    ax1.legend(bbox_to_anchor=(0.5, -0.2), loc='upper center', 
               ncol=min(4, len(df_sorted)), fontsize=9)
    
    # Plot 2: Intron Phase Statistics (Top-right) - STACKED BARS
    ax2 = fig.add_subplot(gs[0, 1])
    
    # Define colors for the stacked bars
    colors_stopless = ['#1f77b4', '#ff7f0e', '#2ca02c']  # Blue, Orange, Green for stopless
    colors_stop_containing = ['#d62728', '#9467bd', '#8c564b']  # Red, Purple, Brown for stop containing
    
    # Set up positions for bars
    y_pos = np.arange(len(df_sorted))
    bar_width = 0.35
    offset = 0.2  # Space between the two bars for each species
    
    # Plot stopless introns (first bar for each species)
    bottom_stopless = np.zeros(len(df_sorted))
    stopless_components = ['stopless_3n', 'stopless_3n1', 'stopless_3n2']
    stopless_labels = ['Phase 0', 'Phase 1', 'Phase 2']
    
    for i, (component, label, color) in enumerate(zip(stopless_components, stopless_labels, colors_stopless)):
        values = df_sorted[component].values
        
        # Plot stacked bars
        bars = ax2.barh(y_pos - offset, values, bar_width, 
                       left=bottom_stopless, color=color, alpha=0.8, 
                       edgecolor='black', label=f'Stopless {label}' if i == 0 else "")
        
        # Add percentage labels in the middle of each segment
        for j, (value, y_pos_val) in enumerate(zip(values, y_pos - offset)):
            if value > 5:  # Only label if segment is large enough
                # Calculate position for label
                x_pos = bottom_stopless[j] + value / 2
                ax2.text(x_pos, y_pos_val, f'{value:.1f}%', 
                        ha='center', va='center', fontsize=8, fontweight='bold',
                        color='white' if value > 20 else 'black')
        
        bottom_stopless += values
    
    # Plot stop-containing introns (second bar for each species)
    bottom_stop_containing = np.zeros(len(df_sorted))
    stop_containing_components = ['stop_containing_3n', 'stop_containing_3n1', 'stop_containing_3n2']
    stop_containing_labels = ['Phase 0', 'Phase 1', 'Phase 2']
    
    for i, (component, label, color) in enumerate(zip(stop_containing_components, stop_containing_labels, colors_stop_containing)):
        values = df_sorted[component].values
        
        # Plot stacked bars
        bars = ax2.barh(y_pos + offset, values, bar_width, 
                       left=bottom_stop_containing, color=color, alpha=0.8, 
                       edgecolor='black', label=f'Stop-Containing {label}' if i == 0 else "")
        
        # Add percentage labels in the middle of each segment
        for j, (value, y_pos_val) in enumerate(zip(values, y_pos + offset)):
            if value > 5:  # Only label if segment is large enough
                # Calculate position for label
                x_pos = bottom_stop_containing[j] + value / 2
                ax2.text(x_pos, y_pos_val, f'{value:.1f}%', 
                        ha='center', va='center', fontsize=8, fontweight='bold',
                        color='white' if value > 20 else 'black')
        
        bottom_stop_containing += values
    
    # REMOVED: x-axis label and grid
    ax2.set_xlabel('')
    ax2.set_title('B: Intron Phase Distribution in Short Introns', fontweight='bold')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(df_sorted['ToLID'])
    ax2.set_xlim(0, 100)
    
    # REMOVED: grid lines
    ax2.grid(False)
    
    # Remove x-axis ticks and tick labels
    ax2.set_xticks([])
    ax2.set_xticklabels([])
    
    # Create custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors_stopless[0], edgecolor='black', alpha=0.8, label='Stopless Phase 0'),
        Patch(facecolor=colors_stopless[1], edgecolor='black', alpha=0.8, label='Stopless Phase 1'),
        Patch(facecolor=colors_stopless[2], edgecolor='black', alpha=0.8, label='Stopless Phase 2'),
        Patch(facecolor=colors_stop_containing[0], edgecolor='black', alpha=0.8, label='Stop-Containing Phase 0'),
        Patch(facecolor=colors_stop_containing[1], edgecolor='black', alpha=0.8, label='Stop-Containing Phase 1'),
        Patch(facecolor=colors_stop_containing[2], edgecolor='black', alpha=0.8, label='Stop-Containing Phase 2')
    ]
    
    ax2.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    
    # Plot 3: Structural Features (Bottom-left) - HORIZONTAL BAR PLOTS
    ax3 = fig.add_subplot(gs[1, 0])
    
    # Updated structural metrics with median values including multiexon transcript median
    structural_metrics = [
        'median_exons_per_transcript', 
        'median_exons_per_multiexon_transcript', 
        'median_cds_length', 
        'median_intron_length'
    ]
    
    # Updated legend labels with normalization info
    metric_labels = [
        'Median Exons per\nTranscript', 
        'Median Exons per\nMultiexon Transcript', 
        'Median CDS Length\n/100 (bp)', 
        'Median Intron Length\n/1000 (bp)'
    ]
    
    colors_structural = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']  # Different colors for each metric
    
    y_pos = np.arange(len(df_sorted))
    height = 0.2  # Reduced height to accommodate 4 metrics
    
    for i, (metric, label, color) in enumerate(zip(structural_metrics, metric_labels, colors_structural)):
        values = df_sorted[metric].values
        
        # Normalize CDS and intron lengths for better visualization
        if metric == 'median_cds_length':
            values = values / 100  # Divide by 100 for scaling
        elif metric == 'median_intron_length':
            values = values / 1000  # Divide by 1000 for scaling
        
        # Create horizontal bars
        bars = ax3.barh(y_pos + i * height, values, height, label=label, 
                       color=color, alpha=0.8, edgecolor='black')
        
        # Add value labels for ALL bars in bold
        for bar in bars:
            width = bar.get_width()
            ax3.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
                    f'{width:.1f}', ha='left', va='center', fontsize=8, fontweight='bold')
    
    ax3.set_xlabel('Value (Normalized)')
    ax3.set_title('C: Structural Features Comparison', fontweight='bold')
    ax3.set_yticks(y_pos + height * 1.5)
    ax3.set_yticklabels(df_sorted['ToLID'])
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax3.grid(True, alpha=0.3, axis='x')
    
    # Plot 4: Transcript Statistics (Bottom-right) - STACKED BARS
    ax4 = fig.add_subplot(gs[1, 1])
    
    # Calculate transcripts with introns
    df_sorted['transcripts_with_introns_count'] = (
        df_sorted['num_transcripts'] - df_sorted['transcripts_without_introns_count']
    )
    
    # Normalize counts by dividing by 1000
    transcripts_with_introns = df_sorted['transcripts_with_introns_count'] / 1000
    transcripts_without_introns = df_sorted['transcripts_without_introns_count'] / 1000
    
    # Plot stacked transcripts bars
    y_pos_counts = np.arange(len(df_sorted))
    height_counts = 0.6  # Single bar height
    
    transcripts_intron_bars = ax4.barh(y_pos_counts, transcripts_with_introns, height_counts,
                                      color='#7f8c8d', alpha=0.8, edgecolor='black', 
                                      label='Transcripts with Introns/1000')
    
    transcripts_no_intron_bars = ax4.barh(y_pos_counts, transcripts_without_introns, height_counts,
                                         left=transcripts_with_introns,
                                         color='#2c3e50', alpha=0.8, edgecolor='black',
                                         label='Transcripts without Introns/1000')
    
    # Add value labels for total transcripts (at the end of stacked bars)
    total_transcripts = (df_sorted['num_transcripts'] / 1000).values
    for i, (y_pos, total) in enumerate(zip(y_pos_counts, total_transcripts)):
        ax4.text(total + 0.5, y_pos, f'{total:.1f}k', 
                ha='left', va='center', fontsize=9, fontweight='bold')
    
    # Add percentage labels for transcripts without introns
    for i, (y_pos, without_introns, total) in enumerate(zip(y_pos_counts, 
                                                           transcripts_without_introns, 
                                                           total_transcripts)):
        if without_introns > 0:  # Only label if there are transcripts without introns
            percentage = (without_introns / total) * 100
            # Position the percentage label in the middle of the "without introns" segment
            left_position = total - without_introns / 2
            ax4.text(left_position, y_pos, f'{percentage:.1f}%', 
                    ha='center', va='center', fontsize=8, color='white', fontweight='bold')
    
    ax4.set_xlabel('Count (thousands)')
    ax4.set_title('D: Transcript Statistics', fontweight='bold')
    ax4.set_yticks(y_pos_counts)
    ax4.set_yticklabels(df_sorted['ToLID'])
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax4.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    return fig

def main():
    # Find all JSON files in current directory
    json_files = glob.glob("*.json")
    
    if not json_files:
        print("No JSON files found in current directory!")
        return
    
    print(f"Found {len(json_files)} JSON files")
    
    # Parse all files
    df = parse_json_files(json_files)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("tab10")
    
    print("Creating comprehensive comparison plot...")
    
    # Create the single comprehensive plot
    fig = create_comprehensive_plot(df)
    
    # Save in multiple formats
    fig.savefig('annotation_comprehensive_comparison.tiff', dpi=300, bbox_inches='tight', format='tiff')
    fig.savefig('annotation_comprehensive_comparison.pdf', bbox_inches='tight', format='pdf')
    fig.savefig('annotation_comprehensive_comparison.png', dpi=300, bbox_inches='tight', format='png')
    
    plt.close(fig)
    
    # Save data as CSV
    df.to_csv('annotation_comparison_data.csv', index=False)
    
    print("Comprehensive comparison plot created successfully!")
    print("Output files:")
    print("  - annotation_comprehensive_comparison.tiff")
    print("  - annotation_comprehensive_comparison.pdf") 
    print("  - annotation_comprehensive_comparison.png")
    print("  - annotation_comparison_data.csv")
    
    # Print verification of parsed ToLIDs
    print("\n=== ToLIDs for Verification ===")
    for _, row in df.sort_values('ToLID').iterrows():
        print(f"ToLID: {row['ToLID']:<20} Species: {row['Species']}")

if __name__ == "__main__":
    main()