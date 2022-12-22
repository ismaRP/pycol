import matplotlib.pyplot as plt
from deamidation.MQreader import EvidenceBatchReader
import numpy as np
import pandas as pd
from pyteomics import mass
from Bio.Align import substitution_matrices
import pycol_functions as pf
import seaborn as sns

main_path = '/home/ismael/palaeoproteomics/'
datapath = main_path + 'MSMSdatasets/'
mspath = datapath + 'mq/'

samples_df = pd.read_csv(datapath + 'all_samples.tsv',
                         index_col=0,
                         sep='\t')
substrates = {'Parchment', 'Hide', 'Skin', 'Leather'}
parch_samples = samples_df[samples_df.Substrate.isin(substrates)]

prots_df = pd.read_csv(datapath + 'proteins.csv',
                       index_col=0,
                       sep=',')

reader = EvidenceBatchReader(
    mspath,
    prot_f=prots_df,
    samples_f=parch_samples,
    fasta_f='',
    include='20201228_IR_ParchOxi',
    ptm=['de', 'hy', ''], aa=['QN', 'P', 'C'], tr='log', int_threshold=150, sf_exp=[])

print('reading...')
mqdata = reader.readBatch()

print('filtering for COL1A1 and COL1A2...')
mqdata.filter_prots(keep_names=['COL1A1', 'COL1A2'])
mqdata.filter_nonassociated_peptides(which='razor')

print('assigning PTMs to proteins...')
mqdata.createPep2Prot()
mqdata.assign_mod2prot(which='razor')

pos_map = mqdata.map_positions(['COL1A1', 'COL1A2'], alnformat='clustal', aamod='PC', return_map=True)

unimod_db = mass.Unimod()
demass = mass.calculate_mass(composition=unimod_db.by_title('Deamidated')['composition'])
hymass = mass.calculate_mass(composition=mass.Composition({'O': 1}))


# for mqrun in mqdata.mqruns:
mqrun = mqdata.mqruns[0]
p_sites = {}
gpp_sites = {}
for sample_name, sample in mqrun.samples.items():
    for prot_id, prot in sample.prot_dict.items():
        seq = mqdata.prot_seqs[prot_id].seq
        pf.calc_hyp_sites(prot, seq, p_sites, gpp_sites)

p_sites = pd.DataFrame.from_dict(
    p_sites, orient='index',
    columns=['int_mod', 'int_unmod', 'position', 'site_type'])
p_sites = p_sites.reset_index()
p_sites = p_sites.assign(frac=lambda x: x.int_mod/(x.int_mod + x.int_unmod))

gpp_sites = pd.DataFrame.from_dict(
    gpp_sites, orient='index',
    columns=['int_mod_second', 'int_unmod_second', 'int_mod_third', 'int_unmod_third',
             'position']
)
gpp_sites = gpp_sites.reset_index()
gpp_sites = gpp_sites.assign(
    frac_second=lambda x: x.int_mod_second/(x.int_mod_second + x.int_unmod_second),
    frac_third=lambda x: x.int_mod_third/(x.int_mod_third + x.int_unmod_third)
)

gpp_sites_long = pd.wide_to_long(
    gpp_sites, stubnames=['int_mod', 'int_unmod', 'frac'], sep='_',
    i=['index', 'position'], j='GPPpos', suffix='\w+')
gpp_sites_long = gpp_sites_long.reset_index()


def count_range(df, r):
    fm = df['frac'].to_numpy()
    counts = [np.sum(fm >= i) for i in r]
    df = pd.DataFrame({'ranges': r, 'counts': counts})
    return df


r = np.arange(0, 1, 0.05)
by_site_type = p_sites.groupby(['site_type'])
p_ranges = by_site_type.apply(count_range, r=r).reset_index()

gpp_byP = gpp_sites_long.groupby('GPPpos')
gpp_ranges = gpp_byP.apply(count_range, r=r).reset_index()


fig, axes = plt.subplots(ncols=2, figsize=(8, 4))
fig.set_tight_layout(True)
sns.scatterplot(x='frac_second', y='frac_third', data=gpp_sites, ax=axes[0], color='black')
axes[0].set_xlabel('GPP second position Hyp fraction')
axes[0].set_ylabel('GPP third position Hyp fraction')

sns.lineplot(x='ranges', y='counts', data=gpp_ranges, hue='GPPpos', ax=axes[1], legend=False)
sns.scatterplot(x='ranges', y='counts', data=gpp_ranges, hue='GPPpos', ax=axes[1])
axes[1].set_xlabel('Hyp fraction')
axes[1].set_ylabel('# positions with >= Hyp fraction')
axes[1].legend(loc='upper right', title='Pro pos')

fig.savefig(main_path + 'MALDI/carla/pycol_results/GPPhydroxylation.pdf')
plt.show()
