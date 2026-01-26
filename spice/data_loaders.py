import os
import re
from pathlib import Path

import pandas as pd
import numpy as np
from tqdm.auto import tqdm

from spice import config, directories
from spice.utils import (CALC_NEW, open_pickle, save_pickle, get_logger, log_debug,
                         calc_telomere_bound_whole_arm_whole_chrom)


logger = get_logger('data_loaders')
CHROMS = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
DATA_LOADERS_DIR = os.path.join(directories['results_dir'], 'data_loaders')

@CALC_NEW()
def load_final_events_df_for_tsg_og(results_dir, skip_assertions=False, remove_multiple_wgds=True, skip_overlaps=False):
    raise NotImplementedError('This function is deprecated; Revisit when I implement TSG/OG analysis again.')

    name = results_dir.split('/')[-1]

    # import here to avoid circular imports
    from spice.pipeline_postprocessing import calc_event_overlaps

    print(f'Loading final events from {os.path.join(results_dir, "final_events.pickle")}')

    final_events_df = open_pickle(os.path.join(results_dir, 'final_events.pickle'))
    print(f'Raw events have: {final_events_df["sample"].nunique()} samples, {final_events_df["id"].nunique()} IDs and {len(final_events_df)} events')

    if remove_multiple_wgds:
        old_len = len(final_events_df)
        with open(f'/Users/tom/phd/projects/signatures/data/{name}/multiple_wgd_samples.txt', 'r') as f:
            multiple_wgd_samples = f.read().splitlines()
        final_events_df = final_events_df.query('sample not in @multiple_wgd_samples').copy().reset_index(drop=True)
        print(f'Removed {old_len - len(final_events_df)} events from {len(multiple_wgd_samples)} samples with multiple WGDs')
        print(f'Events now have: {final_events_df["sample"].nunique()} samples, {final_events_df["id"].nunique()} IDs and {len(final_events_df)} events')

    valid_ids = (final_events_df.groupby('id').size().loc[
        (final_events_df.groupby('id').size() ==
         final_events_df.groupby('id')['events_per_chrom'].first())].index.values)
    if len(valid_ids) < final_events_df['id'].nunique():
        logger.warning(f'Found {final_events_df["id"].nunique() - len(valid_ids)} IDs ({100*(final_events_df["id"].nunique() - len(valid_ids)) / final_events_df["id"].nunique():.4f}%) with inconsistent number of events')
        final_events_df = final_events_df.query('id in @valid_ids').copy()
        print('Removed invalid IDs, where the number of events did not match "events_per_chrom"')
        print(f'Events now have: {final_events_df["sample"].nunique()} samples, {final_events_df["id"].nunique()} IDs and {len(final_events_df)} events')


    final_events_df = final_events_df.join(load_centromeres(extended=True), on='chrom')

    final_events_df[
        ['centromere_bound_l', 'centromere_bound_r', 'telomere_bound_l',
         'telomere_bound_r', 'telomere_bound', 'whole_arm', 'whole_chrom']] = np.stack(calc_telomere_bound_whole_arm_whole_chrom(
             final_events_df, return_left_and_right=True), axis=1)
    
    # Adjust start/end, especially important for chrX where the telomere assignment is off
    final_events_df.loc[final_events_df['telomere_bound_l'], 'start'] = 0
    final_events_df.loc[final_events_df['telomere_bound_r'], 'end'] = final_events_df.loc[final_events_df['telomere_bound_r'], 'chrom_length']

    final_events_df['whole_arm'] = final_events_df.eval('(telomere_bound_l and centromere_bound_r) or (telomere_bound_r and centromere_bound_l)')
    final_events_df.loc[final_events_df.query('whole_chrom').index, 'whole_arm'] = False

    final_events_df['centromere_bound'] = np.logical_or(final_events_df['centromere_bound_l'].values, final_events_df['centromere_bound_r'].values)
    final_events_df['whole_centromere'] = np.logical_and(final_events_df['centromere_bound_l'].values, final_events_df['centromere_bound_r'].values)
    # assert that the number of entries per ID is correct before removing whole centromeres and re-assigning events per chrom
    assert skip_assertions or (final_events_df.groupby('id').size() == final_events_df.groupby('id')['events_per_chrom'].first()).all()
    
    # Check whether any event is within 1Mbp of the centromere and remove them
    # Note that observed centromeres are removed in create_features_pipeline
    print('Remove whole centromere events and events within 1Mbp of the centromere')
    centromeres = load_centromeres(extended=True)
    final_events_df = (final_events_df
        .drop(columns=['centro_start', 'centro_end'], errors='ignore')
        .join(centromeres, on='chrom'))
    centromeres_pad = load_centromeres(extended=True, pad=5e6).rename(
        columns={'centro_start': 'centro_start_pad', 'centro_end': 'centro_end_pad'})
    final_events_df = (final_events_df
        .drop(columns=['centro_start_pad', 'centro_end_pad'], errors='ignore')
        .join(centromeres_pad, on='chrom'))
    final_events_df['inside_centromere'] = final_events_df.eval('start>=centro_start_pad-2 and end<=centro_end_pad+2')    
    
    final_events_df = final_events_df.query('not whole_centromere and not inside_centromere').drop(columns=['whole_centromere', 'inside_centromere']).copy()
    assert len(final_events_df) > 0
    final_events_df = final_events_df.drop(columns='events_per_chrom', errors='ignore').join(final_events_df.groupby('id').size().to_frame('events_per_chrom'), on='id')
    print(f'Events now have: {final_events_df["sample"].nunique()} samples, {final_events_df["id"].nunique()} IDs and {len(final_events_df)} events')

    short_chroms = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']
    final_events_df['short_arm'] = False
    final_events_df.loc[final_events_df.query('chrom in @short_chroms').index, 'short_arm'] = True
    final_events_df.loc[final_events_df.query('chrom in @short_chroms and telomere_bound_r and start <= centro_end_pad+2').index, 'whole_chrom'] = True

    # Events that are within 0.95-1.05 of the arm size are considered whole arm
    telomere_bound_events = final_events_df.query('telomere_bound and not whole_chrom and not whole_arm').copy()
    telomere_bound_events['left_bound'] = telomere_bound_events.eval('start <= 100000')
    telomere_bound_events['right_bound'] = telomere_bound_events.eval('end >= chrom_length - 100000')
    telomere_bound_events.loc[telomere_bound_events['chrom'].isin(short_chroms), 'left_bound'] = telomere_bound_events.loc[telomere_bound_events['chrom'].isin(short_chroms)].eval('start <= centro_end')
    telomere_bound_events['arm_size'] = telomere_bound_events['centro_start']
    telomere_bound_events.loc[telomere_bound_events['right_bound'], 'arm_size'] = telomere_bound_events.loc[telomere_bound_events['right_bound'], 'chrom_length'] - telomere_bound_events.loc[telomere_bound_events['right_bound'], 'centro_end']
    telomere_bound_events.loc[telomere_bound_events['chrom'].isin(short_chroms), 'arm_size'] = telomere_bound_events.loc[telomere_bound_events['chrom'].isin(short_chroms), 'chrom_length'] - telomere_bound_events.loc[telomere_bound_events['chrom'].isin(short_chroms), 'centro_end']
    telomere_bound_events['within_arm'] = telomere_bound_events['width'] < telomere_bound_events['arm_size']
    telomere_bound_events['width_norm'] = telomere_bound_events['width'] / telomere_bound_events['arm_size']
    whole_chrom_indices = telomere_bound_events.query('(left_bound and right_bound)').index.values
    whole_arm_indices = telomere_bound_events.query('not (left_bound and right_bound) and width_norm > 0.95 and width_norm < 1.05').index.values
    final_events_df.loc[whole_chrom_indices, ['whole_chrom']] = True
    final_events_df.loc[whole_arm_indices, ['whole_arm']] = True

    final_events_df['events_per_chrom_bin'] = pd.cut(final_events_df['events_per_chrom'], bins=[0, 1, 5, 10, 20, 30, 40])
    final_events_df = final_events_df.drop(columns='has_whole_chrom', errors='ignore').join(final_events_df.groupby('id')['whole_chrom'].any().to_frame('has_whole_chrom'), on='id')
    final_events_df = final_events_df.drop(columns='has_whole_arm', errors='ignore').join(final_events_df.groupby('id')['whole_arm'].any().to_frame('has_whole_arm'), on='id')

    final_events_df['pos'] = final_events_df.apply(
        lambda x: 'whole_chrom' if x['whole_chrom'] else 'whole_arm' if x['whole_arm'] else
        'centromere_bound' if x['centromere_bound'] else 'telomere_bound' if x['telomere_bound'] else 'internal', axis=1)
    final_events_df['pos_detail'] = final_events_df.apply(
        lambda x: 'whole_chrom' if x['whole_chrom'] else 'whole_arm' if x['whole_arm'] else
        'centromere_bound_l' if x['centromere_bound_l'] else 'centromere_bound_r' if x['centromere_bound_r'] else
        'telomere_bound_l' if x['telomere_bound_l'] else 'telomere_bound_r' if x['telomere_bound_r'] else 'internal', axis=1)

    final_events_df['nr_segments'] = final_events_df['diff'].str.rfind('1') - final_events_df['diff'].str.find('1') + 1
    final_events_df['nr_segments_without_loh'] = final_events_df['diff'].str.count('1')
    assert skip_assertions or (final_events_df['nr_segments'] >= final_events_df['nr_segments_without_loh']).all(), final_events_df.query('nr_segments < nr_segments_without_loh')[['diff', 'nr_segments', 'nr_segments_without_loh']]
    assert skip_assertions or final_events_df.query('nr_segments != nr_segments_without_loh')['loh_shortened_event'].all()
    assert skip_assertions or (final_events_df.groupby('id')['solved'].agg(list).map(lambda x: len(np.unique(x))) == 1).all(), 'Some IDs have multiple solved values'
    assert skip_assertions or (final_events_df.groupby('id').size() == final_events_df.groupby('id')['events_per_chrom'].first()).all()

    final_events_df['chrom_id'] = final_events_df['sample'] + ':' + final_events_df['chrom'].astype(str)
    final_events_df['signed_diff'] = final_events_df['diff']
    final_events_df.loc[final_events_df['type'].isin(['loh', 'pre-loss', 'loss', 'post-loss']), 'signed_diff'] = (
        final_events_df.loc[final_events_df['type'].isin(['loh', 'pre-loss', 'loss', 'post-loss']), 'diff'].map(lambda x: x.replace('1', '-1')))
    final_events_df['type_diff'] = final_events_df['type'] + ':' + final_events_df['diff']

    final_events_df['chrom_int'] = final_events_df['chrom'].apply(lambda x: x.replace('chr', ''))
    final_events_df['chrom_int'] = pd.Categorical(final_events_df['chrom_int'], categories=[str(i) for i in range(1, 23)] + ['X', 'Y'], ordered=True)

    # Filter bad telomere regions on chrX
    final_events_df = final_events_df.query('chrom != "chrX" or not telomere_bound or (end > 28e5 and start < 1548e5)').copy()

    if not skip_overlaps:
        print('Calculating event overlaps')
        final_events_df = final_events_df.reset_index(drop=True)
        event_overlaps = calc_event_overlaps(
            final_events_df,
            calc_new_filename=os.path.join(directories['results_dir'], f'analysis/{name}/events_overlaps.pickle'),
            calc_new_force_new=True)
        final_events_df = pd.concat([final_events_df, event_overlaps], axis=1)

    # important for some downstream analysis that requires unique indices
    final_events_df = final_events_df.reset_index(drop=True)

    assert skip_assertions or not final_events_df.isna().sum().any()

    print(f'Done! Final events have: {final_events_df["sample"].nunique()} samples, {final_events_df["id"].nunique()} IDs and {len(final_events_df)} events')
    return final_events_df


def resolve_data_file(return_raw=False) -> str:
    """Resolve the chromosome segments file path.
    """
    from spice import config, directories
    name = config.get('name')
    data_dir = config['directories']['data_dir']
    orig = config['input_files']['copynumber']
    processed = os.path.join(data_dir, f"{name}_processed.tsv")
    if not return_raw:
        orig = orig.replace('.tsv', '_split.tsv')
        processed = os.path.join(data_dir, f"{name}_processed_split.tsv")

    # Prefer processed split file if available, else original
    cur_file = orig
    if processed and os.path.exists(processed):
        cur_file = processed
    
    if not os.path.isabs(cur_file):
        cur_file = os.path.join(directories['base_dir'], cur_file)
    log_debug(logger, f"Resolved chrom_segments_file: {cur_file}")
    return cur_file


def load_final_events():
    results_dir = os.path.join(directories['results_dir'], config['name'])
    if not os.path.exists(os.path.join(results_dir, 'final_events.tsv')):
        raise FileNotFoundError(f"final_events.tsv not found in dir {results_dir}. Run SPICE event inference first")
    final_events_df = pd.read_csv(
        os.path.join(results_dir, 'final_events.tsv'), sep='\t', dtype={'cn': str, 'diff': str})
    return final_events_df


def load_segmentation(size=None, data_loaders_dir_top=DATA_LOADERS_DIR):
    # import here to avoid circular imports
    from spice.event_analysis.segmentation import create_segmentation
    cur_filename = os.path.join(data_loaders_dir_top, 'segmentations', f'segmentation_{int(size)}.pickle')
    if not os.path.exists(cur_filename):
        logger.warning(f'File not found: {cur_filename} -> creating segmentation with size {size}')
        if size is not None:
            segmentation = create_segmentation(size)
            save_pickle(segmentation, cur_filename)
        else:
            raise ValueError('Segmentation file not found and size is None')
    else:
        segmentation = open_pickle(cur_filename, fail_if_nonexisting=True)

    return segmentation


def load_segmented(name, which='events', seg='100kbp', cancer_type='pancancer', wgd='allwgd', data_loaders_dir_top=DATA_LOADERS_DIR):

    assert wgd in ['allwgd', 'nowgd', 'wgd', 'pre', 'post']

    if which == 'events':
        cur_file = os.path.join(data_loaders_dir_top, f'{name}/segmented_data', f'final_events_in_segmentation_{cancer_type}_{wgd}_{seg}.pickle')
    elif which == 'events_single':
        cur_file = os.path.join(data_loaders_dir_top, f'{name}/segmented_data', f'final_events_in_segmentation_{cancer_type}_{wgd}_{seg}_single.pickle')
    elif which == 'data_full':
        cur_file = os.path.join(data_loaders_dir_top, f'{name}/segmented_data/cn_data_segmentation_{cancer_type}_{wgd}_{seg}.pickle')
    elif which == 'data':
        cur_file = os.path.join(data_loaders_dir_top, f'{name}/segmented_data/cn_data_filtered_segmentation_{cancer_type}_{wgd}_{seg}.pickle')
    else:
        raise ValueError(f'Unknown which {which}')

    cur = open_pickle(cur_file)
    if cur is None:
        logger.warning(f'File not found: {cur_file}. Recreate using the script "analysis_scripts/script_common_segmentation.py"')
    
    return cur


@CALC_NEW()
def load_genes(filter_out_duplicates=True):
    genes = (pd.read_csv(config['input_files']['genes'], sep='\t')
         .query('chrom in @CHROMS')
         .rename({'name': 'name_encode', 'name2': 'name'}, axis=1)
            )
    # Note: txStart refers to the transcription start site (TSS) of a gene, while cdsStart indicates the start position of the coding sequence (CDS)
    genes = pd.concat([
        genes.groupby('name')['chrom'].first().to_frame('chrom'),
        genes.groupby('name')['txStart'].min().to_frame('start'),
        genes.groupby('name')['txEnd'].max().to_frame('end'),
        # genes.groupby('name')['cdsStart'].min().to_frame('start'),
        # genes.groupby('name')['cdsEnd'].max().to_frame('end'),
        genes.groupby('name')['chrom'].nunique().to_frame('n_chroms')],
        axis=1).sort_values(['chrom', 'start'])
    if filter_out_duplicates:
        genes = genes.query('n_chroms == 1').drop(columns='n_chroms')
    genes = genes.reset_index()
    genes['width'] = genes.eval('end - start')
    genes['pos'] = genes.eval('start + width / 2')
    genes = genes.query('(end - start) > 100')
    genes = genes.loc[~genes.reset_index()['name'].str.contains('ENSG00').values]

    chrom_lengths = load_chrom_lengths()

    breakpoint_dict_1Mbp = {
        chrom: np.append((np.arange(0, chrom_lengths.loc[chrom], 1_000_000))[:-1], chrom_lengths.loc[chrom])
        for chrom in CHROMS}
    breakpoint_dict_100kbp = {
        chrom: np.append((np.arange(0, chrom_lengths.loc[chrom], 100_000))[:-1], chrom_lengths.loc[chrom])
        for chrom in CHROMS}

    cur = []
    for chrom in CHROMS:
        cur_genes = genes.query('chrom == @chrom').copy()
        cur_genes['1Mbp_bin'] = pd.cut(cur_genes['pos'], bins=breakpoint_dict_1Mbp[chrom])
        cur_genes['1Mbp_start'] = cur_genes['1Mbp_bin'].apply(lambda x: x.left).astype(int)
        cur_genes['1Mbp_end'] = cur_genes['1Mbp_bin'].apply(lambda x: x.right).astype(int)
        cur_genes['1Mbp_pos'] = 3/2*cur_genes['1Mbp_start'] + 1/2 * cur_genes['1Mbp_end']
        cur_genes['100kbp_bin'] = pd.cut(cur_genes['pos'], bins=breakpoint_dict_100kbp[chrom])
        cur_genes['100kbp_start'] = cur_genes['100kbp_bin'].apply(lambda x: x.left).astype(int)
        cur_genes['100kbp_end'] = cur_genes['100kbp_bin'].apply(lambda x: x.right).astype(int)
        cur_genes['100kbp_pos'] = 3/2*cur_genes['100kbp_start'] + 1/2 * cur_genes['100kbp_end']
        cur.append(cur_genes)
    genes = pd.concat(cur)
    genes = genes.set_index('name')

    return genes


@CALC_NEW()
def calc_summary_from_events_df(events_df, chrom_segments):
    summary = events_df.groupby('id')[['sample', 'chrom', 'allele', 'solved', 'n_paths', 'events_per_chrom', 'chrom_length', 'actual_chrom_length']].first()
    summary = summary.join(((chrom_segments.set_index('id')[['cn']] == 0).groupby('id')).sum().rename({'cn': 'lohs'}, axis=1), on='id')
    summary['sample_chrom_id'] = summary.index.map(lambda x: ':'.join(x.split(':')[:-1]))
    summary = summary.join((chrom_segments.groupby('id')['cn'].max() <= 1).to_frame('is_0_1'), on='id')
    summary = summary.join((summary.groupby('sample_chrom_id').size() == 2).to_frame('has_sister_id'), 
                on='sample_chrom_id')
    summary = summary.join(events_df.groupby('id')['has_wgd'].first().to_frame('WGD'), on='id')
    assert events_df['SV_overlap'].max() <= 1
    summary['SV_overlaps'] = events_df.groupby('id')['SV_overlap'].sum().loc[summary.index]

    return summary


def load_chrom_lengths(assembly='hg19'):
    if assembly != 'hg19':
        raise ValueError('Only hg19 is supported for now')
    cytobands = pd.read_csv(config['input_files']['cytoband'], sep='\t', header=None)
    cytobands.columns = ['chrom', 'start', 'end', 'name', 'gieStain']
    chrom_lengths = cytobands.groupby('chrom')['end'].max().to_frame('chrom_length')

    chrom_lengths = chrom_lengths.reset_index()
    chrom_lengths['chrom'] = format_chromosomes(chrom_lengths['chrom'])
    chrom_lengths = chrom_lengths.set_index('chrom').sort_index()

    return chrom_lengths


def load_cn_tsv(input_file, alleles=['cn_a', 'cn_b']):
    data = pd.read_csv(input_file, sep='\t', dtype={allele: 'int64' for allele in alleles})
    data = (data
            .infer_objects()
            .rename({'chr': 'chrom', 'sample': 'sample_id', 'major_cn': 'cn_a', 'minor_cn': 'cn_b'}, axis=1)
            [['sample_id', 'chrom', 'start', 'end', 'cn_a', 'cn_b']]
            # .set_index(['sample_id', 'chrom', 'start', 'end'])
            )

    for allele in alleles:
        data.loc[data[allele] > 8, allele] = 8

    data['chrom'] = format_chromosomes(data['chrom'])
    # data = data.set_index(['sample_id', 'chrom', 'start', 'end'])
    # data = data.sort_index()

    data['width'] = data.eval('end - start')

    return data


def load_centromeres(extended=True, observed=False, pad=None):
    '''Create file using create_observed_centromeres_and_telomeres'''

    assert not (extended and observed), 'Cannot have both extended and observed centromeres'
    centromeres = pd.read_csv(os.path.join(Path(__file__).parent.parent, 'objects', 'centromeres_ext.tsv' if extended else ('centromeres_observed.tsv' if observed else 'centromeres.tsv')), sep='\t',  header=[0, 1] if observed else [0], index_col=0)

    if pad is not None:
        centromeres['centro_start'] = np.maximum(centromeres['centro_start'] - pad, 0)
        centromeres['centro_end'] = centromeres['centro_end'] + pad
    return centromeres


def load_telomeres_observed():
    '''Create file using create_observed_centromeres_and_telomeres'''
    telomeres_observed = pd.read_csv(os.path.join(Path(__file__).parent.parent, 'objects', 'telomeres_observed.tsv'), sep='\t',  header=[0, 1], index_col=0)

    return telomeres_observed


def create_observed_centromeres_and_telomeres(final_events_df):
    # import here to avoid circular imports
    from spice.features import DEFAULT_SEGMENT_SIZE_DICT, DEFAULT_LENGTH_SCALES
    centromeres = load_centromeres(extended=False)

    actual_centro_pos = pd.DataFrame(index=CHROMS[:-1], columns=pd.MultiIndex.from_product([['small', 'mid1', 'mid2', 'large'], ['centro_start', 'centro_end']]))
    actual_telomere_pos = pd.DataFrame(index=CHROMS[:-1], columns=pd.MultiIndex.from_product([['small', 'mid1', 'mid2', 'large'], ['chrom_start', 'chrom_end']]))
    for cur_chrom in tqdm(CHROMS[:-1]):
        for cur_length_scale in ['small', 'mid1', 'mid2', 'large']:
            cur_features = [f'{pm} int {i}' for i in DEFAULT_LENGTH_SCALES[cur_length_scale] for pm in ['+', '-']]
            centro_center = centromeres.loc[cur_chrom].mean()
            if centromeres.loc[cur_chrom, 'centro_start'] == 0 or cur_chrom in ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']:
                cur_start = 0
            else:
                cur_start = final_events_df.query('pos == "internal" and chrom == @cur_chrom and end < @centro_center and feature in @cur_features')['end'].max()
                if np.isnan(cur_start):
                    cur_start = centromeres.loc[cur_chrom, 'centro_start']
                else:
                    cur_start = int(np.floor(cur_start/DEFAULT_SEGMENT_SIZE_DICT[cur_length_scale])*DEFAULT_SEGMENT_SIZE_DICT[cur_length_scale])

            actual_centro_pos.loc[cur_chrom, (cur_length_scale, 'centro_start')] = cur_start
            cur_end = final_events_df.query('pos == "internal" and chrom == @cur_chrom and start > @centro_center and feature in @cur_features')['start'].min()
            cur_end = int(np.ceil(cur_end/DEFAULT_SEGMENT_SIZE_DICT[cur_length_scale])*DEFAULT_SEGMENT_SIZE_DICT[cur_length_scale])
            actual_centro_pos.loc[cur_chrom, (cur_length_scale, 'centro_end')] = cur_end

            actual_telomere_pos.loc[cur_chrom, (cur_length_scale, 'chrom_start')] = final_events_df.query('pos == "internal" and chrom == @cur_chrom and feature in @cur_features')['start'].min()
            actual_telomere_pos.loc[cur_chrom, (cur_length_scale, 'chrom_end')] = final_events_df.query('pos == "internal" and chrom == @cur_chrom and feature in @cur_features')['end'].max()

    actual_telomere_pos.to_csv('/Users/tom/phd/projects/signatures/spice/objects/telomeres_observed.tsv', sep='\t')
    actual_centro_pos.to_csv('/Users/tom/phd/projects/signatures/spice/objects/centromeres_observed.tsv', sep='\t')


def load_chrom_lengths():
    chrom_lengths = pd.read_csv(os.path.join(directories['base_dir'], config['input_files']['chrom_lengths']), sep='\t').set_index('chrom')['chrom_length']
    return chrom_lengths


def load_arm_lengths():  
    centromeres = load_centromeres(extended=False)
    chrom_lens = load_chrom_lengths()

    arm_lens = []
    for chrom in centromeres.index:
        if chrom in ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']:
            continue
        arm_lens.append([chrom + 'p', centromeres.loc[chrom, 'centro_start']])
        arm_lens.append([chrom + 'q', chrom_lens.loc[chrom] - centromeres.loc[chrom, 'centro_end']])
    arm_lens = pd.DataFrame(arm_lens, columns=['arm', 'arm_length']).set_index('arm')['arm_length']

    return arm_lens


def load_fragile_sites(which='humcfs'):
    '''
    humcfs downloaded from: https://webs.iiitd.edu.in/raghava/humcfs/download.html

    Mainly the same as https://www.sciencedirect.com/science/article/pii/S1874939907001794 which was used in pcawg
    '''
    if which == 'humcfs':
        fragile_sites_dir = '/Users/tom/phd/projects/signatures/data/fragile_sites/humcfs_fragile.gff3'

        fragile_sites = []
        for x in os.listdir(fragile_sites_dir):
            if not x.endswith('.gff3'):
                continue
            try:
                cur = pd.read_csv(os.path.join(fragile_sites_dir, x), sep='\t', comment='#', header=None,
                            names=['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
                fragile_sites.append(cur)
            except:
                pass
        fragile_sites = pd.concat(fragile_sites)
        fragile_sites[['name', 'inducer', 'frequency', 'band']] = fragile_sites['attributes'].str.split(';', expand=True)
        fragile_sites[['name', 'inducer', 'band']] = fragile_sites[['name', 'inducer', 'band']].map(lambda x: x.split('=')[1])
        fragile_sites[['frequency']] = fragile_sites[['frequency']].map(lambda x: x.split('=')[1].lower() if '=' in x else x)
        fragile_sites[['start', 'end']] = fragile_sites[['start', 'end']].astype(int)
        fragile_sites['pos'] = fragile_sites.eval('(start + end) / 2').astype(int)
        fragile_sites['width'] = fragile_sites.eval('end - start').astype(int)
        fragile_sites.loc[fragile_sites['frequency']=='commom', 'frequency'] = 'common'
        fragile_sites = fragile_sites.query('(frequency == "common" or frequency == "rare") and width < 1e8')
        fragile_sites = fragile_sites.drop(columns=['score', 'strand', 'phase', 'attributes', 'source'])
    else:
        raise ValueError('unknown fragile sites provider')

    return fragile_sites


def format_chromosomes(ds):
    '''copied from medicc.tools'''

    ds = ds.astype('str')
    pattern = re.compile(r"(chr|chrom)?(_)?(0)?((\d+)|X|Y)", flags=re.IGNORECASE)
    matches = ds.apply(pattern.match)
    matchable = ~matches.isnull().any()
    if matchable:
        newchr = matches.apply(lambda x:f"chr{x[4].upper():s}")
        numchr = matches.apply(lambda x:int(x[5]) if x[5] is not None else -1)
        chrlevels = np.sort(numchr.unique())
        chrlevels = np.setdiff1d(chrlevels, [-1])
        chrcats = [f"chr{i}" for i in chrlevels]
        if 'chrX' in list(newchr):
            chrcats += ['chrX',]
        if 'chrY' in list(newchr):
            chrcats += ['chrY',]
        newchr = pd.Categorical(newchr, categories=chrcats)
    else:
        logger.warning("Could not match the chromosome labels. Rename the chromosomes according chr1, "
                      "chr2, ... to avoid potential errors."
                      "Current format: {}".format(ds.unique()))
        newchr = pd.Categorical(ds, categories=ds.unique())
    assert not newchr.isna().any(), "Could not reformat chromosome labels. Rename according to chr1, chr2, ..."
    return newchr