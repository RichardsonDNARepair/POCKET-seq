# MOTREC figure 2 creation script
# 4-25-26
# Richardson Lab

import pandas as pd
from upsetplot import UpSet, from_indicators
import matplotlib.pyplot as plt
# ---------------------------
# load data
# ---------------------------

# load pocket-seq
vegfa_motrec = pd.read_csv("VEGFA_2_MOTREC.txt", sep="\t")
socs3_g4_motrec = pd.read_csv("SOCS3_g4_MOTREC.txt", sep="\t")
socs3_g1_motrec = pd.read_csv("SOCS3_g1_MOTREC.txt", sep="\t")
hek293_s4_motrec = pd.read_csv("HEK293_s4_MOTREC.txt", sep="\t")

# Filtering parameters
def filter_motrec(df):
    return df[(df['fold_enrichment'] >= 4) & (df['-log10(pvalue)'] >= 2)]

vegfa_motrec     = filter_motrec(vegfa_motrec)
socs3_g4_motrec  = filter_motrec(socs3_g4_motrec)
socs3_g1_motrec  = filter_motrec(socs3_g1_motrec)
hek293_s4_motrec = filter_motrec(hek293_s4_motrec)

# load CRISPOR (in-vitro prediction)
vegfa_CRISPOR = pd.read_excel("VEGFA_2_CRISPORofftargets_hg19.xls", engine="xlrd", skiprows=8)
socs3_g4_CRISPOR = pd.read_excel("SOCS3g4_CRISPORofftargets_hg19.xls", engine="xlrd", skiprows=8)
socs3_g1_CRISPOR = pd.read_excel("SOCS3g1_CRISPORofftargets_hg19.xls", engine="xlrd", skiprows=8)
hek293_s4_CRISPOR = pd.read_excel("HEK293site4_CRISPORofftargets_hg19.xls", engine="xlrd", skiprows=8)

# load guide-seq
guide_seq = pd.read_excel("guide_seq.xlsx")
vegfa_guide_seq = guide_seq[guide_seq['Targetsite'] == 'VEGFA_site2']
hek293_guide_seq = guide_seq[guide_seq['Targetsite'] == 'HEK293_sgRNA4']

# load discover-seq+ (only Vegfa site 2)
discover_seq = pd.read_excel("dseq+.xlsx", sheet_name="K562_Vs2_KU")
discover_seq[['chr', 'coords']] = discover_seq['# region string'].str.split(':', expand=True)
discover_seq[['start', 'end']] = discover_seq['coords'].str.split('-', expand=True)
discover_seq = discover_seq.drop(columns='coords')
discover_seq['start'] = discover_seq['start'].astype(int)
discover_seq['end'] = discover_seq['end'].astype(int)

# load Kuscu (HA_dCas9 - HEK293site4)
kuscu = pd.read_excel("Kuscu2014.xlsx", sheet_name="sgRNA4")
kuscu.columns = ['OT_id', 'chr', 'start', 'end', 'peak_score', 'mismatches',
    'overlapping_genes', 'strand',
    'pos1','pos2','pos3','pos4','pos5','pos6','pos7','pos8','pos9','pos10',
    'pos11','pos12','pos13','pos14','pos15','pos16','pos17','pos18','pos19','pos20',
    'pam1','pam2','pam3']
kuscu = kuscu.dropna(subset=['chr'])
seq_cols = [f'pos{i}' for i in range(1, 21)]
pam_cols = ['pam1', 'pam2', 'pam3']
kuscu['guide_sequence'] = kuscu[seq_cols].apply(lambda r: ''.join(r.dropna().astype(str)), axis=1)
kuscu['pam'] = kuscu[pam_cols].apply(lambda r: ''.join(r.dropna().astype(str)), axis=1)
kuscu['start'] = kuscu['start'].astype(int)
kuscu['end'] = kuscu['end'].astype(int)
kuscu['peak_score'] = kuscu['peak_score'].astype(float)
kuscu['mismatches'] = kuscu['mismatches'].astype(int)
kuscu = kuscu[['OT_id', 'chr', 'start', 'end', 'peak_score', 'mismatches', 'strand', 'guide_sequence', 'pam', 'overlapping_genes']]

# ---------------------------
# convert to long data format
# ---------------------------

def rank_col(df, score_col, ascending=False):
    return df[score_col].rank(ascending=ascending, method='min').astype(int)

def prep_crispor(df, guide_name):
    df = df.copy()
    df = df.rename(columns={'chrom': 'chr'})
    df['method'] = 'CRISPOR'
    df['guide'] = guide_name
    df['score'] = df['mitOfftargetScore']
    df['rank'] = rank_col(df, 'score', ascending=False)
    df['orientation'] = df['strand']
    df['guide_sequence'] = df['guideSeq']
    return df[['method', 'guide', 'chr', 'start', 'end', 'orientation', 'score', 'rank', 'guide_sequence']]

vegfa_CRISPOR_long     = prep_crispor(vegfa_CRISPOR,     'VEGFA_2')
socs3_g1_CRISPOR_long  = prep_crispor(socs3_g1_CRISPOR,  'SOCS3_g1')
socs3_g4_CRISPOR_long  = prep_crispor(socs3_g4_CRISPOR,  'SOCS3_g4')
hek293_s4_CRISPOR_long = prep_crispor(hek293_s4_CRISPOR, 'HEK293_s4')

def prep_motrec(df, guide_name):
    df = df.copy()
    df['method'] = 'POCKET-seq'
    df['guide'] = guide_name
    df['score'] = df['fold_enrichment']
    df['rank'] = rank_col(df, 'score', ascending=False)
    df['orientation'] = None
    df['guide_sequence'] = None
    return df[['method', 'guide', 'chr', 'start', 'end', 'orientation', 'score', 'rank', 'guide_sequence']]

vegfa_motrec_long     = prep_motrec(vegfa_motrec,     'VEGFA_2')
socs3_g1_motrec_long  = prep_motrec(socs3_g1_motrec,  'SOCS3_g1')
socs3_g4_motrec_long  = prep_motrec(socs3_g4_motrec,  'SOCS3_g4')
hek293_s4_motrec_long = prep_motrec(hek293_s4_motrec, 'HEK293_s4')

def prep_guideseq(df, guide_name):
    df = df.copy()
    df['method'] = 'GUIDE-seq'
    df['guide'] = guide_name
    df = df.rename(columns={
        '#Chromosome': 'chr',
        'Start': 'start',
        'End': 'end',
        'Strand': 'orientation',
        'GUIDE-Seq Reads': 'score',
        'Target_Sequence': 'guide_sequence'
    })
    df['rank'] = rank_col(df, 'score', ascending=False)
    return df[['method', 'guide', 'chr', 'start', 'end', 'orientation', 'score', 'rank', 'guide_sequence']]

vegfa_guideseq_long  = prep_guideseq(vegfa_guide_seq,  'VEGFA_2')
hek293_guideseq_long = prep_guideseq(hek293_guide_seq, 'HEK293_s4')

discover_long = discover_seq.copy()
discover_long['method'] = 'DISCOVER-seq+'
discover_long['guide'] = 'VEGFA_2'
discover_long['score'] = discover_long[' Ctotal']
discover_long['rank'] = rank_col(discover_long, 'score', ascending=False)
discover_long['orientation'] = discover_long[' Cas9 sense']
discover_long['guide_sequence'] = discover_long[' expected target sequence']
discover_long = discover_long[['method', 'guide', 'chr', 'start', 'end', 'orientation', 'score', 'rank', 'guide_sequence']]

kuscu_long = kuscu.copy()
kuscu_long['method'] = 'HA-dCas9'
kuscu_long['guide'] = 'HEK293_s4'
kuscu_long['score'] = kuscu_long['peak_score']
kuscu_long['rank'] = rank_col(kuscu_long, 'score', ascending=False)
kuscu_long['orientation'] = kuscu_long['strand']
kuscu_long = kuscu_long.rename(columns={'guide_sequence': 'guide_sequence'})
kuscu_long = kuscu_long[['method', 'guide', 'chr', 'start', 'end', 'orientation', 'score', 'rank', 'guide_sequence']]

long_table = pd.concat([
    vegfa_CRISPOR_long, socs3_g1_CRISPOR_long, socs3_g4_CRISPOR_long, hek293_s4_CRISPOR_long,
    vegfa_motrec_long, socs3_g1_motrec_long, socs3_g4_motrec_long, hek293_s4_motrec_long,
    vegfa_guideseq_long, hek293_guideseq_long,
    discover_long,
    kuscu_long
], ignore_index=True)

# ---------------------------
# collapse long data to merge
# ---------------------------

d = 50

def build_collapsed_table(long_table, d):

    results = []

    for (guide, chrom), group in long_table.groupby(['guide', 'chr']):
        group = group.sort_values('start').reset_index(drop=True)

        clusters = []
        current_cluster = [0]

        for i in range(1, len(group)):
            prev = group.iloc[current_cluster[-1]]
            curr = group.iloc[i]
            if curr['start'] - d <= prev['end'] + d:
                current_cluster.append(i)
            else:
                clusters.append(current_cluster)
                current_cluster = [i]
        clusters.append(current_cluster)

        for cluster in clusters:
            rows = group.iloc[cluster]
            collapsed = {
                'guide': guide,
                'chr': chrom,
                'start': int(rows['start'].mean()),
                'end': int(rows['end'].mean()),
                'orientation': rows['orientation'].mode()[0] if rows['orientation'].notna().any() else None,
                'Present_CRISPOR': rows['method'].eq('CRISPOR').any(),
                'Present_POCKET': rows['method'].eq('POCKET-seq').any(),
                'Present_GUIDEseq': rows['method'].eq('GUIDE-seq').any(),
                'Present_DISCOVERseq_plus': rows['method'].eq('DISCOVER-seq+').any(),
                'Present_HA_dCas9': rows['method'].eq('HA-dCas9').any(),
                'mean_rank': rows['rank'].mean().round(1),
                'n_methods': rows['method'].nunique(),
            }
            results.append(collapsed)

    return pd.DataFrame(results)

collapsed_table = build_collapsed_table(long_table, d=d)

with pd.ExcelWriter("results/off_target_compare_tables.xlsx", engine="openpyxl") as writer:
    long_table.to_excel(writer, sheet_name="long_table", index=False)
    collapsed_table.to_excel(writer, sheet_name="collapsed_table", index=False)

# ---------------------------
# make upset plots
# ---------------------------

method_cols = ['Present_CRISPOR', 'Present_POCKET', 'Present_GUIDEseq',
               'Present_DISCOVERseq_plus', 'Present_HA_dCas9']

names = {
    'Present_CRISPOR': 'CRISPOR',
    'Present_POCKET': 'POCKET-seq',
    'Present_GUIDEseq': 'GUIDE-seq',
    'Present_DISCOVERseq_plus': 'DISCOVER-seq+',
    'Present_HA_dCas9': 'HA-dCas9'
}

for guide, group in collapsed_table.groupby('guide'):
    active_cols = [c for c in method_cols if group[c].any()]
    subset = group[active_cols].rename(columns=names)
    upset_data = from_indicators(subset)
    fig = plt.figure(figsize=(10, 5))
    upset = UpSet(upset_data, subset_size='count', show_counts=True, sort_by='cardinality')
    upset.plot(fig)
    plt.suptitle(f'Off-target detection: {guide}', fontsize=14, y=1.02)
    plt.savefig(f'results/upset_{guide}.png', dpi=150, bbox_inches='tight')
    plt.close()
