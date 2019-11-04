import pandas as pd
import os
from collections import defaultdict
import scipy.spatial.distance as spdist
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import click


def config_params(font_size=7):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


def snvs_order():
    order = []
    first = ['A', 'C', 'G', 'T']
    pyr = ['C', 'T']
    for p in pyr:
        for mut in first:
            if mut != p:
                for f in first:
                    for f2 in first:
                        comb = '{}[{}>{}]{}'.format(f, p, mut, f2)
                        order.append(comb)
    return order

# read WGS PCAWG signatures
def pcawg_canonical_snvs():
    script_location = os.path.dirname(os.path.realpath(__file__))
    pcawg_sbs_file = '{}/signature_files/sigProfiler_SBS_signatures_2018_03_28.csv'.format(script_location)
    pcawg_snvs = pd.read_csv(pcawg_sbs_file)
    pcawg_snvs.index = pcawg_snvs.apply(lambda x: '{}[{}>{}]{}'.format(
        x['SubType'][0],
        x['SubType'][1],
        x['Type'][-1],
        x['SubType'][2]), axis=1)

    pcawg_snvs = pcawg_snvs.loc[snvs_order()]
    pcawg_snvs.drop(['SubType', 'Type'], axis=1, inplace=True)
    d_pcawg = {}

    for col in pcawg_snvs.columns:
        d_pcawg[col] = pcawg_snvs[col].tolist()

    return d_pcawg, pcawg_snvs.index.tolist()

# read exome PCAWG signatures
def pcawg_canonical_snvs_exomes():
    script_location = os.path.dirname(os.path.realpath(__file__))
    pcawg_sbs_file = '{}/signature_files/signatures.exome.cosmic.v3.may2019.tsv'.format(script_location)
    pcawg_snvs = pd.read_csv(pcawg_sbs_file, sep ='\t')
    pcawg_snvs = pcawg_snvs.T
    pcawg_snvs = pcawg_snvs.loc[snvs_order()]

    d_pcawg = {}
    for col in pcawg_snvs.columns:
        d_pcawg[col] = pcawg_snvs[col].tolist()

    return d_pcawg, pcawg_snvs.index.tolist()

def plot_cosine_similarity(cos_df, ttype, outpath):
    config_params(5.5)
    fig, ax = plt.subplots(1, 1, figsize=(32, len(cos_df.T) * 0.9))

    sns.heatmap(cos_df.T, annot=True, fmt='g', cmap='YlGnBu', ax = ax)
    ax.set_ylim(-0.5, len(cos_df.T) +0.5)

    os.makedirs(os.path.join(outpath, 'processes', ttype), exist_ok=True)
    plt.savefig('{}/processes/{}/{}.cosine_similarity.png'.format(outpath, ttype, ttype), dpi=300,
                bbox_inches='tight')
    plt.savefig('{}/processes/{}/{}.cosine_similarity.svg'.format(outpath, ttype, ttype, ))
    plt.close()


def get_new_names_signatures(similar, signature_similarity_cutoff):
    new_sig = 1
    dict_equi = {}
    for i, row in similar.iterrows():
        if row[1] > signature_similarity_cutoff:
            dict_equi[i] = '{}_{}-{}'.format(i, row[0], round(float(row[1]), 2))
        else:
            dict_equi[i] = '{}_NA'.format(i)
        new_sig += 1

    return dict_equi


def get_similarities_signatures(df_processes, outpath, signature_similarity_cutoff, ttype, exome):
    if exome is True:
        d_pcawg, index = pcawg_canonical_snvs_exomes()
    else:
        d_pcawg, index = pcawg_canonical_snvs()

    # find similar signatures to those reported in SigProfiler
    cos_sim = defaultdict(dict)
    for ix, col in enumerate(df_processes.columns):
        vec1 = df_processes[col].tolist()
        for s, vec2 in d_pcawg.items():
            c = 1 - round(spdist.cosine(vec1, vec2), 3)
            cos_sim[col][s] = c

    # select those with higher similarity per each of the signatures
    cos_df = pd.DataFrame(cos_sim)
    plot_cosine_similarity(cos_df, ttype, outpath)

    index_max = cos_df.idxmax()
    vals = cos_df.max()
    similar = pd.DataFrame(list(zip(index_max, vals)))

    similar.index = cos_df.columns
    similar.to_csv('{}/processes/{}/{}.cosine_similarity_maximum.tsv'.format(outpath, ttype, ttype), sep ='\t',
                   index = True, header = True)

    dict_equi = get_new_names_signatures(similar, signature_similarity_cutoff)

    df_processes_cols = df_processes.columns

    df_processes.columns = [dict_equi[c] for c in df_processes_cols]

    return df_processes


# split into even chunks
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


# function to plot the SNV processes
def plot_snvs(sig, title, outpath, ttype):
    config_params(3)
    fig, axs = plt.subplots(
        nrows=2, ncols=1, figsize=(3.2, 1), gridspec_kw={'height_ratios': [1, 9]}
    )
    order_plot = snvs_order()

    vals = []
    colors = []
    colors_mut = [
        '#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'
    ]
    bot = -0.5
    for ix, c in enumerate(chunks(sig, 16)):
        colors.extend([colors_mut[ix] for _ in c])
        axs[0].barh(1, 16, left=bot, color=colors_mut[ix])
        bot += 16
        vals.extend(c)

    axs[0].set_xlim(-1, 96)
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['bottom'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)

    axs[0].get_yaxis().set_visible(False)
    axs[0].get_xaxis().set_visible(False)

    x = [i for i in range(len(vals))]

    axs[1].bar(x, vals, color=colors, width=0.8, linewidth=0, align='center')
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(
        ['{}{}{}'.format(a[0], a[2], a[-1]) for a in order_plot],
        verticalalignment="center", ha='center', rotation=90, fontsize=2,
        color='grey'
    )
    axs[1].set_ylabel('Relative Probability')

    plt.tight_layout()
    plt.xlim(-1, 96)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)

    plt.setp([axs[1].get_xticklines(), axs[1].get_yticklines()], color='grey')

    axs[1].xaxis.set_ticks_position('none')
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[1].spines[axis].set_linewidth(0.2)

    axs[1].xaxis.set_tick_params(pad=0.5)
    axs[1].yaxis.set_tick_params(pad=0.5, width=0.5)

    prev_pos = 6
    for count_lab, lab in enumerate(['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']):
        color = 'black'
        if count_lab == 1:
            color = 'white'
        axs[0].text(prev_pos, 0.85, lab, color=color, weight='bold')
        prev_pos = prev_pos + 16

    plt.tick_params(axis='both', which='both', bottom=False, left=False)
    plt.suptitle(title, y=1.05)

    os.makedirs(os.path.join(outpath, 'processes', ttype), exist_ok=True)
    plt.savefig('{}/processes/{}/{}.{}.png'.format(outpath, ttype, ttype, title), dpi=300, bbox_inches='tight')
    plt.savefig('{}/processes/{}/{}.{}.svg'.format(outpath, ttype, ttype, title))

    plt.close()


def process_exposures(path_results, K, outpath, columns, ttype):
    extracted_exposures = os.path.join(path_results, 'exposures_fitting_{}'.format(K))
    exp_df = pd.read_csv(extracted_exposures, sep='\t')

    exp_df.index = columns
    os.makedirs(os.path.join(outpath, 'exposures', ttype), exist_ok=True)

    exp_df.to_csv('{}/exposures/{}/{}.exposures.tsv'.format(outpath, ttype, ttype),
                  sep='\t', index=True, header=True)

    norm = exp_df / exp_df.sum()

    config_params(5.5)
    heatmap = sns.clustermap(norm.fillna(0), cmap='YlGnBu', figsize=(len(norm), len(norm.T)*0.08)
                           )
    plt.setp(heatmap.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)

    heatmap.savefig('{}/exposures/{}/{}.heatmap.png'.format(outpath, ttype, ttype), dpi=300)
    heatmap.savefig('{}/exposures/{}/{}.heatmap.svg'.format(outpath, ttype, ttype))

    exp_df.to_csv('{}/exposures/{}/{}.exposures.tsv'.format(outpath, ttype, ttype),
                        sep ='\t', index = True, header = True)
    plt.close()


def process_signatures(path_results, K, outpath, signature_similarity_cutoff, exome):

    tumor_name = os.path.basename(os.path.normpath(path_results))
    extracted_processes = os.path.join(path_results, 'processes_{}'.format(K))
    extracted_stabilities = os.path.join(path_results, 'processesStabAvg_{}'.format(K))

    df_processes = pd.read_csv(extracted_processes, sep='\t')
    df_stabilities = pd.read_csv(extracted_stabilities, sep='\t')

    dic_stability = df_stabilities.loc[0].to_dict()
    new_cols = ['{}_{}'.format(c.split('x')[1], round(dic_stability[c], 2)) for c in df_processes.columns]
    df_processes.columns = new_cols

    df_processes = get_similarities_signatures(df_processes, outpath, signature_similarity_cutoff, tumor_name, exome)
    df_processes.to_csv('{}/processes/{}/{}.processes.tsv'.format(outpath, tumor_name, tumor_name),
                        sep ='\t', index = False, header = True)

    # plot signatures
    for s, sig in df_processes.iteritems():
        plot_snvs(list(sig), s, outpath, tumor_name)

    process_exposures(path_results, K, outpath, df_processes.columns, tumor_name)


@click.command()
@click.option('--path_results',
              type=click.Path(exists=True),
              help="Path that was used as an outpath in SigProfilerJulia extraction",
              required=True)
@click.option('--sigs_active',
              type=click.INT,
              help="Number of signatures active in the samples",
              required=True)
@click.option('--outpath',
              type=click.Path(),
              help="Path that will be used to save the results",
              required=True)
@click.option('--signature_similarity_cutoff',
              type=click.FLOAT,
              default=0.85,
              help="Cutoff to decide whether the signature resembles a previous signature",
              required=True)
@click.option('--exome',
              is_flag = True,
              default=False,
              help="Compare to COSMIC exome-based signatures",
              )
def run(path_results, sigs_active, outpath, signature_similarity_cutoff, exome):
    process_signatures(path_results, sigs_active, outpath, signature_similarity_cutoff, exome)

if __name__ == '__main__':
    run()
