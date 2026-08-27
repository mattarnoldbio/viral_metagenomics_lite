#!/bin/bash -l

blast_results=$1
which_db=$2
rankedlineage=$3
fullnamelineage=$4

set -euo pipefail

rankedlineage_clean="rankedlineage_clean.tsv";
fullnamelineage_clean="fullnamelineage_clean.tsv";

# full diamond blast table with taxonomic information
joined_table="diamond_blastx_results_with_lineage_all_contigs.tsv";

# If you allowed multiple hits per contig, then this file will store the best hit per contig
best_hit="diamond_best_hit_per_contig.tsv";

# Filtered blast table to include only assigned viruses and excludes phages, retroviruses, and picobirnaviruses
filtered="${which_db}_all_virus_hits.csv";

# ---------------------------------------------------------------------------
# Six recognized virus realms (post-March-2025 NCBI taxonomy), plus "Viruses"
# itself as a fallback in case any rows are still populated the old way.
# Field 10 of rankedlineage.dmp is now "domain/realm" -- for cellular
# organisms it holds Bacteria/Archaea/Eukaryota, for viruses it holds one
# of these realm names instead of the literal word "Viruses".
# ---------------------------------------------------------------------------
VIRUS_DOMAINS="Viruses|Riboviria|Duplodnaviria|Monodnaviria|Varidnaviria|Ribozyviria|Adnaviria";

# ---------------------------------------------------------------------------
# 1: Parse rankedlineage.dmp -- fields 2-10 map directly, in order, onto
#    name/species/genus/family/order/class/phylum/kingdom/domain
# ---------------------------------------------------------------------------
awk -F'\t\\|\t' 'BEGIN{OFS="\t"}
{
    gsub(/\t\|$/, "", $10);
    for (i = 2; i <= 10; i++) if ($i == "") $i = "NA";
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
    # taxid, name, species, genus, family, order, class, phylum, kingdom, domain
}' "$rankedlineage" > "$rankedlineage_clean";

# ---------------------------------------------------------------------------
# 2: Parse fullnamelineage.dmp -- ancestor chain, then append the taxon's
#    own name so it matches the R output's format (lineage ends in self)
# ---------------------------------------------------------------------------
awk -F'\t\\|\t' 'BEGIN{OFS="\t"}
{
    name = $2;
    lineage = $3;
    gsub(/\t\|$/, "", lineage);
    gsub(/; $/, "", lineage);
    print $1, lineage "; " name
}' "$fullnamelineage" > "$fullnamelineage_clean";

# ---------------------------------------------------------------------------
# 3: Split qseqid into Sample_id / query_id, extract first taxid, join
#    against both lookup tables, emit the full 29-column table
# ---------------------------------------------------------------------------
gawk -F'\t' 'BEGIN{OFS="\t"}
ARGIND==1 {
    rl_name[$1]=$2; rl_species[$1]=$3; rl_genus[$1]=$4; rl_family[$1]=$5;
    rl_order[$1]=$6; rl_class[$1]=$7; rl_phylum[$1]=$8; rl_kingdom[$1]=$9;
    rl_domain[$1]=$10;
    next
}
ARGIND==2 { fl_lineage[$1]=$2; next }
ARGIND==3 {
    split($1, id, "|");
    sample_id = id[1]; query_id = id[2];

    split($16, taxids, ";");
    t = taxids[1];

    name    = (t in rl_name)    ? rl_name[t]    : "NA";
    species = (t in rl_species) ? rl_species[t] : "NA";
    genus   = (t in rl_genus)   ? rl_genus[t]   : "NA";
    family  = (t in rl_family)  ? rl_family[t]  : "NA";
    order_  = (t in rl_order)   ? rl_order[t]   : "NA";
    class_  = (t in rl_class)   ? rl_class[t]   : "NA";
    phylum  = (t in rl_phylum)  ? rl_phylum[t]  : "NA";
    kingdom = (t in rl_kingdom) ? rl_kingdom[t] : "NA";
    domain  = (t in rl_domain)  ? rl_domain[t]  : "NA";
    lineage = (t in fl_lineage) ? fl_lineage[t] : "NA";

    print sample_id, query_id, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12,
        $13, $14, $15, t, $17, $18,
        name, species, genus, family, order_, class_, phylum, kingdom, domain, lineage
}' "$rankedlineage_clean" "$fullnamelineage_clean" "$blast_results" > "${joined_table}.body";

echo -e "Sample_id\tquery_id\ttarget_id\t%ident\tlength\tevalue\tbitscore\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tqseq_translated\tqlen\tslen\ttaxid\tstitle\tsscinames\tname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain\tlineage" > "$joined_table";
cat "${joined_table}.body" >> "$joined_table";
rm "${joined_table}.body";

echo "Lineage-joined table written: $joined_table";
wc -l "$joined_table";

# ---------------------------------------------------------------------------
# 4: Keep the single best hit per contig (highest bitscore)
# ---------------------------------------------------------------------------
tail -n +2 "$joined_table" \
    | sort -t$'\t' -k1,1 -k2,2 -k7,7gr \
    | awk -F'\t' '!seen[$1 SUBSEP $2]++' > "${best_hit}.body";

head -1 "$joined_table" > "$best_hit";
cat "${best_hit}.body" >> "$best_hit";
rm "${best_hit}.body";

echo "Best-hit-per-contig table written: $best_hit";
wc -l "$best_hit";

# ---------------------------------------------------------------------------
# 5: Filter -- exclude ONLY unassigned viruses and bacteriophage taxa.
#    Domain check now uses VIRUS_DOMAINS (realm names) instead of a
#    single "Viruses" string match.
# ---------------------------------------------------------------------------
EXCLUDE_TAXA="Caudoviricetes|Microviricetes|Efunaviria|Volvereviria|Vidaverviricetes|Ainoaviricetes|Prepoliviricotina|Tectiviridae|Plasmaviridae|Matsushitaviridae|Obscuriviridae|Vinavirales|Leviviricetes|Cystoviridae|Retroviridae|Picobirnaviridae";

tail -n +2 "$best_hit" | awk -F'\t' -v exclude="$EXCLUDE_TAXA" -v virus_domains="$VIRUS_DOMAINS" 'BEGIN{OFS="\t"}
{
    order = $24; domain = $28; lineage = $29; name = $20;

    if (domain !~ ("^(" virus_domains ")$") && lineage !~ /Viruses/) next;
    if (order == "NA" && lineage !~ /viridae/) next;               # unassigned/unclassified at the order level
    # if (lineage ~ exclude) next;                                     # known phage/archaeal-virus taxa
    # if (tolower(name) ~ /phage/ || tolower(lineage) ~ /phage/) next; # keyword fallback

    print
}' > "${filtered}.body";

head -1 "$joined_table" > "$filtered";
cat "${filtered}.body" >> "$filtered";
#rm "${filtered}.body";

echo "Filtered (phage + unassigned removed) hits written: $filtered";
wc -l "$filtered";

# ---------------------------------------------------------------------------
# Clean up intermediary lookup files -- not needed after the join is done
# ---------------------------------------------------------------------------
#rm -f "$rankedlineage_clean" "$fullnamelineage_clean" "$best_hit";