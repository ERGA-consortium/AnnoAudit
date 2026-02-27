import statistics
import argparse
import datetime
import json
import csv
import re

def parse_statistics_output(stats_file):
    basic_stats = {}
    with open(stats_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            key = row[0].strip()
            value = row[1].strip()
            basic_stats[key] = value
    return basic_stats

def parse_busco_output(busco_file, type):
    with open(busco_file, 'r') as f:
        data = json.load(f)
    results = data["results"]
    s = results["one_line_summary"]
    cleaned = s.replace('%', '').replace(f"n:{results['n_markers']}", '').replace(',,', ',')
    if cleaned.endswith(','):
        cleaned = cleaned[:-1]

    if type == "Assem":
        busco_oneline = {"lineage_dataset": data["lineage_dataset"]["name"]}
        busco_oneline["n_markers"] = str(results["n_markers"])
        busco_oneline["BUSCO " + type] = cleaned
    else:
        busco_oneline = {"BUSCO " + type: cleaned}
    busco_complete = str(results["Complete percentage"])

    return busco_oneline, busco_complete

def parse_omark_output(omark_output):
    omark_stats = {}
    with open(omark_output, 'r') as f:
        for line in f:
            # Completeness
            if line.startswith('The clade used was:'):
                omark_stats['OMA_clade'] = line.split(':')[1].strip()
            elif line.startswith('Number of conserved HOGs:'):
                omark_stats['num_conserved_hogs'] = int(line.split(':')[1].strip())
            elif line.startswith('Single:'):
                omark_stats['single'] = line.split(':')[1].strip()
            elif line.startswith('Duplicated:'):
                omark_stats['duplicated'] = line.split(':')[1].strip()
            elif line.startswith('Duplicated, Unexpected:'):
                omark_stats['duplicated_unexpected'] = line.split(':')[1].strip()
            elif line.startswith('Duplicated, Expected:'):
                omark_stats['duplicated_expected'] = line.split(':')[1].strip()
            elif line.startswith('Missing:'):
                omark_stats['missing'] = line.split(':')[1].strip()

            # Consistency
            elif line.startswith('Number of proteins in the whole proteome:'):
                omark_stats['num_proteins_in_proteome'] = line.split(':')[1].strip()
            elif line.startswith('Total Consistent:'):
                omark_stats['total_consistent'] = line.split(':')[1].strip()
            elif line.startswith('Consistent, partial hits:'):
                omark_stats['consistent_partial_hits'] = line.split(':')[1].strip()
            elif line.startswith('Consistent, fragmented:'):
                omark_stats['consistent_fragmented'] = line.split(':')[1].strip()
            elif line.startswith('Total Inconsistent:'):
                omark_stats['total_inconsistent'] = line.split(':')[1].strip()
            elif line.startswith('Inconsistent, partial hits:'):
                omark_stats['inconsistent_partial_hits'] = line.split(':')[1].strip()
            elif line.startswith('Inconsistent, fragmented:'):
                omark_stats['inconsistent_fragmented'] = line.split(':')[1].strip()
            elif line.startswith('Total Contaminants:'):
                omark_stats['total_contaminants'] = line.split(':')[1].strip()
            elif line.startswith('Contaminants, partial hits:'):
                omark_stats['contaminants_partial_hits'] = line.split(':')[1].strip()
            elif line.startswith('Contaminants, fragmented:'):
                omark_stats['contaminants_fragmented'] = line.split(':')[1].strip()
            elif line.startswith('Total Unknown:'):
                omark_stats['total_unknown'] = line.split(':')[1].strip()
    return omark_stats

def parse_brh_output(brh_output):
    brh_stats = {}
    num_lines = 0
    num_splitting_genes_075 = 0
    num_fusion_genes_125 = 0

    with open(brh_output, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            num_lines += 1
            qlen = int(row[2])
            slen = int(row[3])
            if slen < qlen * 0.75:
                num_splitting_genes_075 += 1
            elif slen > qlen * 1.25:
                num_fusion_genes_125 += 1

    brh_stats['num_best_reciprocal_hits'] = str(num_lines)
    brh_stats['num_split_gene_length_<0.75x'] = str(num_splitting_genes_075) + " (" + str(round(num_splitting_genes_075/num_lines*100, 2)) + "%)"
    brh_stats['num_fusion_gene_length_>1.25x'] = str(num_fusion_genes_125) + " (" + str(round(num_fusion_genes_125/num_lines*100, 2)) + "%)"
    return brh_stats

def parse_psauron_output(csv_file_path):
    psauron_score = None
    true_count = 0
    false_count = 0
    scores = []
    
    with open(csv_file_path, newline='') as file:
        reader = csv.reader(file)
        
        next(reader, None) # Skip first line
        score_line = next(reader, None)
        if score_line:
            score_match = re.search(r'psauron score: ([\d\.]+)', score_line[0])
            psauron_score = float(score_match.group(1)) if score_match else None
        next(reader, None) # Skip the third line (header)
        
        # Process remaining lines
        for row in reader:
            if len(row) >= 3:
                if row[1] == "True":
                    true_count += 1
                elif row[1] == "False":
                    false_count += 1
                
                try:
                    scores.append(float(row[2]))
                except ValueError:
                    pass
    
    return {
        "psauron_score": psauron_score,
        "true_count": true_count,
        "false_count": false_count,
        "median_score": statistics.median(scores) if scores else None,
        "max_score": max(scores) if scores else None,
        "min_score": min(scores) if scores else None
    }

def parse_flagstat(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    passed = data["QC-passed reads"]

    return {
        "mapping_rate": str(passed["mapped %"]) + "%",
        "primary_mapping_rate": str(passed["primary mapped %"]) + "%",
        "properly_paired": str(passed["properly paired %"]) + "%"
    }

def parse_feature(file_path):
    feature_stats = {}
    with open(file_path, "r") as f_in:
        for line in f_in:
            data = line.rstrip().split(sep="\t")
            feature_stats[data[0]] = str(data[1])
    return feature_stats

def parse_json(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data

def combine_results(statistics, busco, omark, brh, prot_distribution, transcriptome, featurecount, psauron, intron_info, taxon_info, warnings):
    combined_stats = {**statistics, **busco, **omark, **brh, **prot_distribution, **transcriptome, **featurecount, **psauron, **intron_info, **taxon_info}
    protein_combined = {**brh, **prot_distribution}
    rnaseq_combined = {**transcriptome, **featurecount, **intron_info}

    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %Z")
    

    with open("Evaluation_output.txt", "w") as f_out:
        
        max_key_length = max(len(str(key)) for key in combined_stats)
        max_value_length = max(len(str(value)) for value in combined_stats.values())
        format_string = f"|{{:<{max_key_length}}} | {{:<{max_value_length}}} |"

        print(format_string.format("Taxon Info", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in taxon_info.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("General Statistics", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in statistics.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("BUSCO", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in busco.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("OMARK", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in omark.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("PSAURON", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in psauron.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("Best Reciprocal Hits", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in protein_combined.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("RNASeq", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in rnaseq_combined.items():
            print(format_string.format(key, value), file=f_out)

    summary = {"num_genes": str(statistics["num_genes"]),
               "num_transcripts": str(statistics["num_transcripts"]),
               "num_transcripts_without_introns": str(statistics["num_transcripts_without_introns"]),
               "mean_exons_per_transcript": str(statistics["mean_exons_per_transcript"]),
               "mean_exons_per_multiexon_transcript": str(statistics["mean_exons_per_multiexon_transcript"]),
               "mean_cds_length": str(statistics["mean_cds_length"]),
               "median_cds_length": str(statistics["median_cds_length"]),
               "mean_intron_length": str(statistics["mean_intron_length"]),
               "median_intron_length": str(statistics["median_intron_length"]),
               "BUSCO_annotation": str(busco["BUSCO Annot"]),
               "BUSCO_assembly": str(busco["BUSCO Assem"]),
               "PSAURON_score": str(psauron["psauron_score"]),
               "num_intron_supported": str(rnaseq_combined["num_exact_intron_boundary"])
               }

    combined_stats = {"Taxon Info": taxon_info,
                      "Warnings": warnings,
                      "Summary": summary,
                      "General Statistics": statistics,
                      "BUSCO": busco,
                      "OMArk": omark, 
                      "PSAURON": psauron,
                      "Best Reciprocal Hits": protein_combined, 
                      "RNASeq": rnaseq_combined}
    
    pretty_json_obj = json.dumps(combined_stats, indent=4)
    with open("Evaluation_output.json", "w") as json_out:
        json_out.write(pretty_json_obj)

def main():
    parser = argparse.ArgumentParser(description="Combine and analyze input files.")
    parser.add_argument("-s", "--statistics_output", help="Path to the general statistics file")
    parser.add_argument("-b", "--busco_annot_output", help="Path to the busco annotation output file")
    parser.add_argument("-l", "--busco_assem_output", help="Path to the busco assembly output file")
    parser.add_argument("-o", "--omark_output", help="Path to the omark output file")
    parser.add_argument("-r", "--brh_output", help="Path to the brh output file")
    parser.add_argument("-i", "--idxstats_output", help="Path to the idxstats output file")
    parser.add_argument("-c", "--compare_distribution_output", help="Path to the compare distribution output")
    parser.add_argument("-f", "--feature_output", help="Path to the featureCounts output")
    parser.add_argument("-p", "--psauron_output", help="Path to the output of Psauron")
    parser.add_argument("-k", "--intron_output", help="Path to the output of intron analysis")
    parser.add_argument("-t", "--taxon_info", help="Path to the taxon information")
    parser.add_argument("-w", "--reference_warning", help="Threshold for warning comparison")
    parser.add_argument("--tolid", help="TOLID of the genome")
    args = parser.parse_args()

    statistics_out = parse_statistics_output(args.statistics_output)
    busco_assembly, busco_assem_complete = parse_busco_output(args.busco_assem_output, "Assem")
    busco_annotation, busco_annot_complete = parse_busco_output(args.busco_annot_output, "Annot")
    omark_out = parse_omark_output(args.omark_output)
    brh_out = parse_brh_output(args.brh_output)
    transcriptome_out = parse_flagstat(args.idxstats_output)
    prot_distribution_out = parse_statistics_output(args.compare_distribution_output)
    feature_out = parse_feature(args.feature_output)
    psauron_out = parse_psauron_output(args.psauron_output)
    intron_out = parse_feature(args.intron_output)
    taxon_info = parse_feature(args.taxon_info)
    ref_warning = parse_json(args.reference_warning)

    taxon_info["ToLID"] = args.tolid

    delta_busco_value = round(float(busco_annot_complete) - float(busco_assem_complete),2)
    delta_busco = {"Delta BUSCO": str(delta_busco_value) + "%"}
    busco = {**busco_assembly, **busco_annotation, **delta_busco}

    num_supported_intron = int(intron_out["num_exact_intron_boundary"])
    percent_intron_supported = round(num_supported_intron/int(statistics_out["num_introns"])*100, 2)
    intron_out["num_exact_intron_boundary"] = "{} ({}%)".format(str(num_supported_intron), str(percent_intron_supported))

    warnings = {"1. Number of genes": "No warnings",
               "2. Predicted CDS": "No warnings",
               "3. BUSCO/OMArk completeness": "No warnings",
               "4. Intron length distribution": "No warnings",
               "5. RNA-Seq support": "No warnings"}
    
    def extract_percentage(text):
        pattern = r'\(\s*(\d+\.\d+)\s*%\s*\)|(\d+\.\d+)\s*%'
        match = re.search(pattern, text)
        if match:
            return match.group(1) or match.group(2)
        return None
    
    warning_1 = []
    if int(statistics_out["num_genes"]) > ref_warning["num_genes"]:
        warning_1.append("Number of genes > {} ({})".format(str(ref_warning["num_genes"]), str(statistics_out["num_genes"])))
    percent_transcript_wo_intron = extract_percentage(str(statistics_out["num_transcripts_without_introns"]))
    if float(percent_transcript_wo_intron) > ref_warning["num_transcripts_without_introns"]:
        warning_1.append("Number monoexon transcript > {}% ({}%)".format(str(ref_warning["num_transcripts_without_introns"]), str(percent_transcript_wo_intron)))

    warning_2 = []
    if float(statistics_out["median_cds_length"]) < ref_warning["median_cds_length"]:
        warning_2.append("Median CDS length < {} ({})".format(str(ref_warning["median_cds_length"]), str(statistics_out["median_cds_length"])))
    if float(psauron_out["psauron_score"]) < ref_warning["psauron_score"]:
        warning_2.append("PSAURON score < {} ({})".format(str( ref_warning["psauron_score"]), str(psauron_out["psauron_score"])))
    percent_split_gene = extract_percentage(brh_out["num_split_gene_length_<0.75x"])
    if float(percent_split_gene) > ref_warning["num_split_gene_length_<0.75x"]:
        warning_2.append("Number possible split gene > {}% ({}%)".format(str(ref_warning["num_split_gene_length_<0.75x"]), str(percent_split_gene)))
    percent_fusion_gene = extract_percentage(brh_out["num_fusion_gene_length_>1.25x"])
    if float(percent_fusion_gene) > ref_warning["num_fusion_gene_length_>1.25x"]:
        warning_2.append("Number possible fusion gene > {}% ({}%)".format(str(ref_warning["num_fusion_gene_length_>1.25x"]), str(percent_fusion_gene)))

    warning_3 = []
    if float(busco_annot_complete) < ref_warning["busco_annot"]:
        warning_3.append("BUSCO Annotation completeness < {}% ({}%)".format(str(ref_warning["busco_annot"]), str(busco_annot_complete)))
    if abs(float(delta_busco_value)) > ref_warning["delta_busco"]:
        warning_3.append("Delta BUSCO > {}% ({}%)".format(str(ref_warning["delta_busco"]), str(delta_busco_value)))
    percent_omark_contaminant = extract_percentage(omark_out["total_contaminants"])
    if float(percent_omark_contaminant) > ref_warning["total_contaminants"]:
        warning_3.append("OMArk contaminants > {}% ({}%)".format(str(ref_warning["total_contaminants"]), str(percent_omark_contaminant)))

    warning_4 = []
    percent_short_3n_intron = extract_percentage(statistics_out["short_intron_3n/short_intron"])
    if float(percent_short_3n_intron) > ref_warning["short_intron_3n"]:
        warning_4.append("Excess short 3n intron > {}% ({}%)".format(str(ref_warning["short_intron_3n"]), str(percent_short_3n_intron)))
    percent_stopless_short_3n_intron = extract_percentage(statistics_out["stopless_short_intron_3n/stopless_short_intron"])
    if float(percent_stopless_short_3n_intron) > ref_warning["stopless_short_intron_3n"]:
        warning_4.append("No deficit of stopless short intron 3n, the percentage is > {}% ({}%)".format(str(ref_warning["stopless_short_intron_3n"]), str(percent_stopless_short_3n_intron)))

    warning_5 = []
    mapping_rate = extract_percentage(transcriptome_out["mapping_rate"])
    if float(mapping_rate) < ref_warning["mapping_rate"]:
        warning_5.append("Genome mapping rate < {}% ({}%)".format(str(ref_warning["mapping_rate"]), str(mapping_rate)))
    percent_gene_sup = extract_percentage(feature_out["num_gene_supported"])
    if float(percent_gene_sup) < ref_warning["num_gene_supported"]:
        warning_5.append("Percent gene supported < {}% ({}%)".format(str(ref_warning["num_gene_supported"]), str(percent_gene_sup)))
    percent_exon_sup = extract_percentage(feature_out["num_exon_supported"])
    if float(percent_exon_sup) < ref_warning["num_exon_supported"]:
        warning_5.append("Percent exon supported < {}% ({}%)".format(str(ref_warning["num_exon_supported"]), str(percent_exon_sup)))
    if float(percent_intron_supported) < ref_warning["num_intron_supported"]:
        warning_5.append("Percent intron supported < {}% ({}%)".format(str(ref_warning["num_intron_supported"]), str(percent_intron_supported)))
    
    if len(warning_1) > 0:
        warnings["1. Number of genes"] = warning_1
    if len(warning_2) > 0:
        warnings["2. Predicted CDS"] = warning_2
    if len(warning_3) > 0:
        warnings["3. BUSCO/OMArk completeness"] = warning_3
    if len(warning_4) > 0:
        warnings["4. Intron length distribution"] = warning_4
    if len(warning_5) > 0:
        warnings["5. RNA-Seq support"] = warning_5

    combine_results(statistics_out, busco, omark_out, brh_out, prot_distribution_out, transcriptome_out, feature_out, psauron_out, intron_out, taxon_info, warnings)

if __name__ == '__main__':
    main()