#!/usr/bin/env python3

import pandas as pd
import subprocess
import argparse
import glob

def add_unique_id(df):
    df['unique_id'] = (df['sample'] + "-" + df['contig']).astype('string')
    return df

def load_hit_data(file_path, filestring, reads=False):
    if not reads:
        headers = ["sample", "contig", "species", "genus", "family", "score", "multi", "length"]
    else:
        headers = ["sample", "contig", "species", "genus", "family", "score", "length"]

    diamond_virus_hits = pd.read_csv(file_path + "/diamond" + filestring + "_all_virus_hits.csv", names=headers)
    blast_virus_hits = pd.read_csv(file_path + "/blastn" + filestring + "_all_virus_hits.csv", names=headers)
    blast_non_virus_hits = pd.read_csv(file_path + "/blastn" + filestring + "_non_virus_hits.csv", names=headers)
    [diamond_virus_hits, blast_virus_hits, blast_non_virus_hits] = [add_unique_id(df) for df in [diamond_virus_hits, blast_virus_hits, blast_non_virus_hits]]
    return diamond_virus_hits, blast_virus_hits, blast_non_virus_hits

vertebrate_families = ['anelloviridae', 'circoviridae', 'caliciviridae', 'reoviridae', 'matonaviridae', 'spinareoviridae', 'tobaniviridae', 'hepadnaviridae', 'nyamiviridae', 'hepeviridae', 'rhabdoviridae', 'sedoreoviridae', 'papillomaviridae', 'pneumoviridae', 'birnaviridae', 'astroviridae', 'vilyaviridae', 'poxviridae', 'virgavirida', 'nodaviridae', 'hantaviridae', 'phenuiviridae', 'filoviridae', 'paramyxoviridae', 'adintoviridae', 'kolmioviridae', 'naryaviridae', 'flaviviridae', 'arenaviridae', 'redondoviridae', 'orthomyxoviridae', 'nairoviridae', 'asfarviridae', 'adenoviridae', 'herpesviridae', 'togaviridae', 'ViPR', 'peribunyaviridae', 'smacoviridae', 'bornaviridae', 'coronaviridae', 'parvoviridae', 'orthoherpesviridae', 'arteriviridae', 'polyomaviridae', 'picornaviridae', 'none']

non_vert_species = "|".join(['bacteriophage', 'phage', 'caudoviricetes','myoviridae','siphoviridae','levi-like','lenarviricota','pahexavirus', 'mycovirus', 'mycoalphavirus', 'plateaulakevirus', 'dolichocephalovirinae'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter Diamond hits')
    parser.add_argument('-d','--data_dir', type=str, help='Path to the directory containing the hit files', default=None)
    parser.add_argument('--diamond_cutoff', type=float, help='Diamond cutoff score', default=0.0)
    parser.add_argument('--blast_cutoff', type=float, help='Blast cutoff score', default=0.0)
    parser.add_argument('--min_contig_length', type=int, help='Minimum contig length', default=250)
    parser.add_argument('--is_reads', action='store_true')
    parser.set_defaults(is_reads=False)
    parser.add_argument('--filter_vertebrate', action='store_true')
    parser.set_defaults(filter_vertebrate=False)
    args = parser.parse_args()

    is_reads = args.is_reads
    if not is_reads:
        print("Processing results of diamond on contigs")
        hits = "/all_contigs_hits.fa"
        filestring = ""
    elif is_reads:
        print("Processing results of diamond on unmapped reads")
        hits = "/all_reads_hits.fa"
        filestring = "_read_level"
    

    diamond_virus_hits, blast_virus_hits, blast_non_virus_hits = load_hit_data(args.data_dir, filestring, reads=is_reads)    

    diamond_virus_hits = diamond_virus_hits[diamond_virus_hits['score'] <= args.diamond_cutoff]
    diamond_virus_hits = diamond_virus_hits[diamond_virus_hits['length'] >= args.min_contig_length]
    blast_non_virus_hits = blast_non_virus_hits[blast_non_virus_hits['score'] <= args.blast_cutoff]
    blast_non_virus_hits = blast_non_virus_hits[blast_non_virus_hits['length'] >= args.min_contig_length]

    blast_filtered = set(diamond_virus_hits.unique_id.to_list()) - set(blast_non_virus_hits.unique_id.to_list())
    blast_removed = set(diamond_virus_hits.unique_id.to_list()) - set(blast_filtered)

    print("num non viral blastn hits: ", len(set(blast_non_virus_hits.unique_id.to_list()))) 
    print("diamond hits removed due to blastn results: ", len(blast_removed))

    blast_removed_diamond_df = diamond_virus_hits[diamond_virus_hits['unique_id'].isin(blast_removed)]
    blast_filtered_diamond_df = diamond_virus_hits[diamond_virus_hits['unique_id'].isin(blast_filtered)]

    removed_hits = pd.merge(blast_removed_diamond_df, blast_non_virus_hits, on='unique_id', how='outer', sort=False, suffixes=('_diamond', '_blast'))

    if args.min_contig_length != 250:
        filestring += "_min_length_" + str(args.min_contig_length)

    removed_hits.drop("unique_id", axis=1).to_csv(args.data_dir + "/diamond" + filestring + "_false_positives.csv", index=False)

    blast_filtered_diamond_df = blast_filtered_diamond_df.drop("unique_id", axis=1)

    contigs = [line for line in open(args.data_dir + hits)]

    sequences = []

    for hit in blast_filtered_diamond_df.iterrows():
        hit = hit[1]
        contig_id = '>' + hit['sample'] + '|' + hit.contig
        # print("contig_id: ", contig_id)
        # print("Data dir:", args.data_dir)
        # print("Hits:", hits)
        contig = subprocess.run(['grep', contig_id, args.data_dir + hits, '-A', '1'], capture_output=True, text=True).stdout.split('\n')[1]
        sequences.append(contig)

    blast_filtered_diamond_df['sequence'] = sequences


    if args.filter_vertebrate:
        print("Filtering out vertebrate hits")
        non_vertebrate_viral_hits = blast_filtered_diamond_df[~blast_filtered_diamond_df['family'].str.lower().isin(vertebrate_families)]
        non_vertebrate_viral_hits = pd.concat([non_vertebrate_viral_hits,  blast_filtered_diamond_df[blast_filtered_diamond_df['species'].str.lower().str.contains(non_vert_species)]], ignore_index=True)
        blast_filtered_diamond_df = blast_filtered_diamond_df[blast_filtered_diamond_df['family'].str.lower().isin(vertebrate_families)]
        blast_filtered_diamond_df = blast_filtered_diamond_df[~blast_filtered_diamond_df['species'].str.lower().str.contains(non_vert_species)]
        non_vertebrate_viral_hits.drop_duplicates().to_csv(args.data_dir + "/final" + filestring + "non_vertebrate_viral_hits.csv", index=False)
        blast_filtered_diamond_df.drop_duplicates().to_csv(args.data_dir + "/final" + filestring + "vertebrate_viral_hits.csv", index=False)
    else:
        blast_filtered_diamond_df.to_csv(args.data_dir + "/final" + filestring + "_viral_hits.csv", index=False)
        