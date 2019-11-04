from glob import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import click
import matplotlib.patches as mpatches
import os

def config_params(font_size=7):
    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = font_size
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'


def get_summary_results(path_results):

    all_stab_avg = glob(os.path.join(path_results,'avgStability_*'))
    all_stab_specific = glob(os.path.join(path_results, 'processesStabAvg_*'))
    all_err_avg = glob(os.path.join(path_results, 'avgReconstructionErrorPercentage_*'))

    results = defaultdict(dict)
    for f in all_stab_avg:
        number_signatures = int(f.split('_')[-1])
        with open(f, 'rt') as infile:
            for line in infile:
                results['stab_avg'][number_signatures] = float(line.rstrip())

    for f in all_err_avg:
        number_signatures = int(f.split('_')[-1])
        with open(f, 'rt') as infile:
            for line in infile:
                results['stab_err'][number_signatures] = float(line.rstrip())

    for f in all_stab_specific:
        number_signatures = int(f.split('_')[-1])
        df = pd.read_csv(f, sep='\t')
        results['avg_specific'][number_signatures] = list(df.loc[0].values)

    return results


def plot_signatures(results, path_results):

    config_params(5)
    fig, ax = plt.subplots(1, 1, figsize=(3, 1.5))

    ax2 = ax.twinx()
    sorted_numbers = sorted(results['stab_avg'].keys())
    lines_avg = []
    lines_err = []
    xticks = []

    for signature in sorted_numbers:

        ax.scatter(signature, results['stab_avg'][signature], color='darkred', s = 8)
        ax.scatter([signature for _ in range(signature)], results['avg_specific'][signature], color='darkred', alpha=0.2, s = 8)
        ax2.scatter(signature, results['stab_err'][signature], color='orange', s = 8)
        lines_avg.append(results['stab_avg'][signature])
        lines_err.append(results['stab_err'][signature])
        xticks.append(signature)

    ax.set_ylabel('Stability')
    ax2.set_ylabel('Reconstruction Error')
    ax.plot(xticks, lines_avg, color='darkred')
    ax2.plot(xticks, lines_err, color='orange')
    ax.set_xlabel('Signatures Active')

    pop_a = mpatches.Patch(color='darkred', label='Average Signatures Stability')
    pop_b = mpatches.Patch(color='orange', label='Average Reconstruction Error')

    plt.legend(handles=[pop_a, pop_b], bbox_to_anchor=(1.1, 1.25))
    plt.savefig('{}/signatures_analysis.png'.format(path_results), dpi=600,  bbox_inches='tight')
    plt.close()


@click.command()
@click.option('--path_results',
              type=click.Path(exists=True),
              help="Input data",
              required=True)
def run(path_results):
    results = get_summary_results(path_results)
    plot_signatures(results, path_results)

if __name__ == '__main__':
    run()
