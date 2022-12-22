from multiprocessing import Pool
import sys
from itertools import islice
import numpy as np
import pandas as pd
import re
from pyteomics import parser
from pyteomics import mass
from collections import deque
import itertools as it
from functools import partial
import pycol_tree as pt


def parallelize_job(ldata, f, n_cores, more_args=None):
    pool = Pool(n_cores)
    if more_args is not None:
        f = partial(f, **more_args)
    pl = pool.map(f, ldata)
    pool.close()
    pool.join()
    return pl


# --------------------------------------------------------------------
# Functions for theoretical markers generation

def cleave(sequence, rule, missed_cleavages=0, min_length=None, semi=False):
    """
    returns tuple (sequence, 1-based position, # misscleavages)
    """
    peptides = []
    ml = missed_cleavages + 2
    n_miss_cl = list(range(missed_cleavages, -1, -1))
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    if min_length is None:
        min_length = 1
    cl = 1
    for i in it.chain([x.end() for x in rule.finditer(sequence)]):
        cleavage_sites.append(i)
        if cl < ml:
            cl += 1
        for j in trange[:cl - 1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq and len(seq) >= min_length:
                peptides.append((i - len(seq) + 1, seq, n_miss_cl[j]))
                if semi:
                    for k in range(min_length, len(seq) - 1):
                        peptides.append((i - len(seq[:k]) + 1, seq[:k], n_miss_cl[j]))
                    for k in range(1, len(seq) - min_length + 1):
                        peptides.append((i - len(seq[k:]) + 1, seq[k:], n_miss_cl[j]))
    return peptides


def modseqs_uniq_comb(pept, var_ptm, varpos_ptm, fixed_ptm, fixedpos_ptm):
    pass
    # ptms are {'label1': ['X', 'Y', ...]} or {'label1': [POS1, POS2, ...]}
    # Go through labels and count how many on pept
    # Then generate all combinations


def generate_modseqs(pepts, mass_lim, chain_lim, var_ptm, uniq_comb=True):
    n_ptm = {}  # {'hy': 3, 'de': 2}
    aa_set = set()  # {'Q', 'N', 'P'}
    for ptm, aa_list in var_ptm.items():
        n_ptm[ptm] = 0
        aa_set.update(aa_list)
    massmin, massmax = mass_lim
    if chain_lim is not None:
        start, end = chain_lim
    else:
        start, end = 0, 1e6
    pept_df = []
    for pos, p, mc in pepts:
        if pos < start or (pos + len(p) - 1) > end:
            continue
        pos = pos - start + 1
        ptm_combs = set()
        peptmass = mass.fast_mass(p)
        if not (massmin <= peptmass <= massmax):
            continue
        # Calculate all mod combinations
        modseqs = parser.isoforms(
            p, variable_mods=var_ptm
        )
        for modseq in modseqs:
            p_modseq = parser.parse(modseq, split=True)
            # Re-initializa n_ptm
            for k in n_ptm:
                n_ptm[k] = 0
            for i, aa in enumerate(p_modseq):
                if len(aa) > 1:
                    ptm = aa[0]
                    aa = aa[1]
                    n_ptm[ptm] += 1
            add_modseq = True
            if uniq_comb:
                comb = ''
                for ptm in sorted(n_ptm.keys()):
                    comb += str(n_ptm[ptm])
                if comb in ptm_combs:
                    add_modseq = False
                else:
                    ptm_combs.add(comb)

            if add_modseq:
                modpept_mass = mass.fast_mass2(modseq, ion_type='M+H')
                modpept_mass = round(modpept_mass, 4)
                new = [p, pos, modpept_mass, mc]
                for ptm in sorted(n_ptm.keys()):
                    new.append(n_ptm[ptm])
                pept_df.append(new)
    return pept_df


def generate_msms_modseqs(pepts, mass_lim, chain_lim, var_ptm, pos_map, pos_mod, seqname,
                          mod_threshold=(0.1, 0.9), uniq_comb=True):
    n_ptm = {}  # {'hy': 3, 'de': 2}
    pept_n_ptm = {}
    aa_set = set()  # {'Q', 'N', 'P'}
    for ptm, aa_list in var_ptm.items():
        n_ptm[ptm] = 0
        pept_n_ptm[ptm] = []
        aa_set.update(aa_list)
    massmin, massmax = mass_lim
    if chain_lim is not None:
        start, end = chain_lim
    else:
        start, end = 0, 1e6
    pept_df = []
    tot_pepts = 0
    for pos, p, mc in pepts:
        if pos < start or (pos + len(p) - 1) > end:
            continue
        peptmass = mass.fast_mass(p)
        if not (massmin <= peptmass <= massmax):
            continue
        # Calculate all mod combinations
        modseqs = parser.isoforms(
            p, variable_mods=var_ptm
        )
        # Re-initialize pept_n_ptm
        for k in pept_n_ptm:
            pept_n_ptm[k] = []
        for modseq in modseqs:
            tot_pepts += 1
            p_modseq = parser.parse(modseq, split=True)
            # Re-initializa n_ptm
            for k in n_ptm:
                n_ptm[k] = 0
            for i, aa in enumerate(p_modseq):
                if len(aa) == 1:
                    ptm = None
                    aa = aa[0]
                else:
                    ptm = aa[0]
                    aa = aa[1]
                    n_ptm[ptm] += 1
                    if n_ptm[ptm] in pept_n_ptm[ptm] and uniq_comb:
                        break
                    else:
                        pept_n_ptm[ptm].append(n_ptm[ptm])
                if aa not in aa_set:  # aa is not investigated
                    continue
                global_pos = pos + 1 + i
                mapped = pos_map[seqname].get(global_pos, 0)
                avg_mod = pos_mod.get(aa + str(mapped), mod_threshold[0])
                if avg_mod > mod_threshold[1]:
                    # This mod is forced
                    if ptm is None:
                        break
                elif mod_threshold[0] <= avg_mod <= mod_threshold[1]:
                    pass  # Both mod or not mod allowed
                elif avg_mod < mod_threshold[0]:
                    if ptm is not None:
                        break
            else:
                modpept_mass = mass.fast_mass2(modseq, ion_type='M+H')
                modpept_mass = round(modpept_mass, 3)
                new = [p, pos, modpept_mass, mc]
                for ptm in sorted(n_ptm.keys()):
                    new.append(n_ptm[ptm])
                pept_df.append(new)
    return pept_df, tot_pepts


def generate_peptides(rec, seq_annot=None, osregex=None, pos_map=None, pos_mod=None, var_ptm=None, missed_cleavages=1,
                      mod_threshold=(0.1, 0.9), mass_lim=(800, 3500), uniq_comb=True, accession_sep='|', verbose=False):
    if seq_annot is None and osregex is None:
        sys.exit('Please specify either a data frame with sequence annotations or \
                  a regular expression to extract them from the fast header')
    if var_ptm is None:
        var_ptm = {'hyd': ['P'], 'deam': ['Q', 'N']}
    n_ptm = {}  # {'hy': 3, 'de': 2}
    pept_n_ptm = {}
    aa_set = set()  # {'Q', 'N', 'P'}
    for ptm, aa_list in var_ptm.items():
        n_ptm[ptm] = 0
        pept_n_ptm[ptm] = []
        aa_set.update(aa_list)

    # Get org name and taxid
    acc = rec.name
    if accession_sep is not None:
        acc = acc.split(accession_sep)[1]
    sequence = str(rec.seq)
    if seq_annot is not None:
        org = seq_annot.loc[acc, 'species']
        taxid = seq_annot.loc[acc, 'taxid']
        name = seq_annot.loc[acc, 'protein_name']
        chain_lim = (seq_annot.loc[acc, 'start'], seq_annot.loc[acc, 'end'])
    else:
        m = osregex.search(rec.description)
        org = m.group(1)
        taxid = int(m.group(2))
        chain_lim = None
        name = None

    enz_re = parser.expasy_rules['trypsin']
    enz_re = re.compile(enz_re)
    pepts = cleave(sequence, enz_re, missed_cleavages=missed_cleavages)
    # pepts = parser.icleave(sequence, enz_re, missed_cleavages=1)
    if pos_map is not None and pos_mod is not None:
        pept_df, tot_pepts = generate_msms_modseqs(
            pepts, mass_lim, chain_lim, var_ptm, pos_map, pos_mod,
            acc, mod_threshold, uniq_comb)
    else:
        pept_df = generate_modseqs(pepts, mass_lim, chain_lim, var_ptm, uniq_comb)

    columns = ['seq', 'seqpos', 'mass1', 'missed.cleaves']
    for ptm in sorted(n_ptm.keys()):
        columns.append('n' + ptm)
    pept_df = pd.DataFrame(
        pept_df,
        columns=columns
    )
    pept_df['org'] = org
    pept_df['taxid'] = taxid
    pept_df['name'] = name
    pept_df['accession'] = acc
    pept_df.drop_duplicates(inplace=True)

    if verbose and pos_map is not None and pos_mod is not None:
        print('Sequence: {}'.format(acc))
        msg = "Found {} peptides in total, "
        msg += "of which {} are MS2 confirmed.\n"
        print(msg.format(tot_pepts, pept_df.shape[0]))
    if verbose and pos_map is None and pos_mod is None:
        print('Sequence: {}'.format(acc))
        msg = "Found {} peptides in total.\n"
        print(msg.format(pept_df.shape[0]))
    return pept_df


def group_mass(df, field='mass1'):
    bin_n = 1
    groups = []
    for gr in df.groupby(field):
        df = gr[1]
        df['bin_n'] = bin_n
        bin_n += 1
        groups.append((gr[0], df))
    return groups


def group_bin(df, size=0.3):
    min_mass = df['mass1'].min()
    max_mass = df['mass1'].max()
    masses = df['mass1'].to_numpy()

    seq = np.arange(min_mass, max_mass, size)
    binidx = np.digitize(masses, seq)
    df['bin_n'] = binidx

    return df.groupby('bin_n')


def group_cluster(df, maxdist=0.3):
    df = df.sort_values('mass1')
    current_bin = 0
    bincol = []
    last_mass = 0
    for row in df.itertuples():
        m = row.mass1
        if abs(m - last_mass) > maxdist:
            current_bin = current_bin + 1
        bincol.append(current_bin)
        last_mass = m
    df['bin_n'] = bincol
    return df.groupby('bin_n')


def agg_pepts(gr, tree, taxid_list):
    df = gr[1]
    seq, ndeam, nhyd = gr[0]
    markers_taxid_list = df['taxid'].unique()
    lca = pt.get_lca(tree, markers_taxid_list)
    org = tree.nodes[lca]['name']

    if df.shape[0] > 1:
        agg = True
    else:
        agg = False
    row1 = df.iloc[:1, :].copy()
    row1['taxid'] = lca
    row1['org'] = org
    row1['aggregated'] = agg
    for t in taxid_list:
        row1[str(t)] = False
    for _, t in df['taxid'].iteritems():
        row1[str(t)] = True
    return pd.DataFrame(row1)


def flag_multiseq(gr):
    df = gr[1]
    seqs = df['seq']
    if seqs.min() != seqs.max():
        df['multiseq'] = True
    else:
        df['multiseq'] = False
    return df


def calc_hyp_sites(prot, seq, p_sites, gpp_sites):
    for ptm in prot.ptms.values():
        pos = ptm[2]  # Used to retrieve neighbor aas
        mapped_pos = ptm[3]  # Used to id the mod
        if pos < 2 or pos >= len(seq) or not (ptm[0] == 'C' or ptm[0] == 'P'):
            continue
        s12 = seq[pos - 3] + seq[pos - 2]
        s2 = seq[pos - 2]
        s3 = seq[pos]
        p_mapped_id = ptm[0] + str(mapped_pos)
        if s12 == 'GP':  # GPP third
            gpp_mapped_pos = mapped_pos - 2
            p_pos = 'GPP-third'
            gpp_mapped_id = 'GPP' + str(gpp_mapped_pos)
        elif s2 == 'G' and s3 == 'P':  # GPP second
            gpp_mapped_pos = mapped_pos - 1
            p_pos = 'GPP-second'
            gpp_mapped_id = 'GPP' + str(gpp_mapped_pos)
        else:
            p_pos = s12 + 'P' + s3
            gpp_mapped_id = ''
            gpp_mapped_pos = 0

        if p_mapped_id not in p_sites:
            p_sites[p_mapped_id] = [0, 0, mapped_pos, p_pos]
        p_sites[p_mapped_id][0] += ptm[5]
        p_sites[p_mapped_id][1] += ptm[4]

        if gpp_mapped_id not in gpp_sites:
            gpp_sites[gpp_mapped_id] = [0, 0, 0, 0, gpp_mapped_pos]
        if p_pos == 'GPP-second':
            gpp_sites[gpp_mapped_id][0] += ptm[5]
            gpp_sites[gpp_mapped_id][1] += ptm[4]
        elif p_pos == 'GPP-third':
            gpp_sites[gpp_mapped_id][2] += ptm[5]
            gpp_sites[gpp_mapped_id][3] += ptm[4]
    gpp_sites.pop('', None)


def calc_nmod(pept):
    nhyd, ndeam = 0, 0
    seq = pept.sequence.replace('GPC', 'GPP')
    for ptm in pept.ptms:
        if ptm[1] == 'C' or ptm[1] == 'P(hy)':
            nhyd += 1
        elif ptm[1] == 'Q(de)' or ptm[1] == 'N(de)':
            ndeam += 1
    return nhyd, ndeam, seq


def find_ms2pepts(row, seq_dic, ptms=None):
    if ptms is None:
        ptms = ['ndeam', 'nhyd']
    seq = row['seq']
    if row['multiseq'] is True or seq not in seq_dic:
        row['pept_ms2conf'] = False
        for ptm in ptms:
            row[ptm + '_ms2conf'] = False
        return row
    else:
        row['pept_ms2conf'] = True
        for ptm in ptms:
            if sum(row[ptm] == seq_dic[seq][ptm]) >= 1:
                row[ptm+'_ms2conf'] = True
            else:
                row[ptm + '_ms2conf'] = False
        return row
