"""Extra preprocessing module.

This module was previously a standalone script. It now exposes a
callable `main(...)` function so it can be invoked from the SPICE CLI.
The general logic and behavior are preserved.
"""

# Imports
import os

import pandas as pd
import numpy as np
import medicc
import fstlib
from medicc.core import create_standard_fsa_dict_from_data

from spice import config, directories, data_loaders
from spice.logging import log_debug, get_logger
from spice.preprocessing.preprocessing import (
    fill_gaps_cnsistent_wrapper, merge_neighbours_mod, infer_wgd_status,
    fill_telomere_nans, get_breaks_mod, main_aggregate_quiet
)
from spice.preprocessing.preprocessing import get_or_infer_wgd_status, get_or_infer_xy_status
from spice.event_inference.fst_assets import nowgd_fst, get_diploid_fsa, T_forced_WGD, SYMBOL_TABLE

logger = get_logger('preprocessing')


def main(
    unique_chroms: bool = False,
    total_cn: bool = False,
    skip_phasing: bool = False,
    skip_centromeres: bool = False,
):
    """Run the extra preprocessing pipeline.

    Parameters
    - name: Project/dataset name used for output paths
    - unique_chroms: Only keep unique chromosomes
    - total_cn: Treat input as total copy-number (single cn column)
    - skip_phasing: Skip MEDICC2 phasing
    - skip_centromeres: Skip centromere binning/unification
    - wgd: WGD inference method, one of {'pcawg', 'spectrum'}
    """
    name = config['name']

    cn_columns = ['total_cn'] if total_cn else ['cn_a', 'cn_b']
    cn_columns_phased = ['cn_a'] if total_cn else ['cn_a', 'cn_b']
    diploid_fsa = get_diploid_fsa(total_copy_numbers=total_cn)

    log_debug(logger, f'{"Total" if total_cn else "Haplotype-specific"} copynumber data with CN columns: {cn_columns}')
    log_debug(logger, f'Saving processed data to {directories["data_dir"]}')

    chrom_lens_df = data_loaders.load_chrom_lengths()

    data = data_loaders.load_raw_copy_number_data(config['input_files']['copynumber'])[['sample_id', 'chrom', 'start', 'end', 'cn_a', 'cn_b']]

    # assert data.index.get_level_values('sample_id').nunique() == len(drews_whitelist), (data.index.get_level_values('sample_id').nunique(), len(drews_whitelist))
    data.loc[data['start'] == 1, 'start'] = 0 # change such that start begins at 0
    data[['start', 'end']] = data[['start', 'end']].astype(int)
    log_debug(logger, f"Found {data['sample_id'].nunique()} samples with a total of {len(data.drop_duplicates(['sample_id', 'chrom']))} unique chromosomes and {len(data)} entries")
    # data = data.set_index(['sample_id', 'chrom', 'start', 'end'])
    all_sample_chrom_ids = data.groupby(['sample_id', 'chrom'], observed=True).size().index.sort_values()

    wgd_status = get_or_infer_wgd_status(data=data, total_cn=total_cn)
    nowgd_samples = wgd_status.loc[~wgd_status.values].index
    wgd_samples = wgd_status.loc[wgd_status.values].index
    
    # Resolve XY (male) samples from file or by inference
    xy_status = get_or_infer_xy_status(data)
    xy_samples = xy_status.loc[xy_status.values].index

    ## Remove minor CN of chrX and chrY for male samples (should only affect a single sample)
    if not total_cn and xy_samples is not None:
        cur_index = data.eval('(sample_id in @xy_samples) and (chrom == "chrX" or chrom == "chrY")').values
        cur = (data.loc[cur_index][['cn_a', 'cn_b']].min(axis=1) > 0).sum()
        log_debug(logger, f"Removing {cur} segments, where male patients have minor CN values for X or Y chroms")
        data.loc[cur_index, 'cn_a'] = data.loc[cur_index, ['cn_a', 'cn_b']].sum(axis=1)
        data.loc[cur_index, 'cn_b'] = 0

    # Cap copynumber at 8 and filter small segments
    if total_cn:
        log_debug(logger, f'Capping copy numbers at max value 8: {data.eval("total_cn > 8").sum()} entries are affected')
        data.loc[data['total_cn'] > 8, 'total_cn'] = 8
    else:
        log_debug(logger, f'Capping copy numbers at max value 8: {data.eval("cn_a > 8").sum()} entries are affected')
        data.loc[data['cn_a'] > 8, 'cn_a'] = 8
        data.loc[data['cn_b'] > 8, 'cn_b'] = 8

    log_debug(logger, f'Removing small segments (1kb): {data.eval(("(end-start) < 1e3")).sum()} entries are affected')
    data = data.query('(end-start) > 1e3').copy().reset_index(drop=True)
    total_cn_sum = data['total_cn'].sum() if total_cn else data[['cn_a', 'cn_b']].sum().sum()

    # set end to chrom length if it is larger
    data = data.join(chrom_lens_df, on='chrom')
    data['end'] = data[['end', 'chrom_length']].min(axis=1)
    data = data.drop(columns='chrom_length')
    assert data.groupby(['sample_id', 'chrom', 'end']).size().max() == 1, 'There were multiple bins with ending larger than chromosome length'

    # Replace gaps with NaN
    gaps = ((data.iloc[1:]['start'].values - data.iloc[:-1]['end'].values)
                * (data.iloc[:-1]['sample_id'].values == data.iloc[1:]['sample_id'].values))
    log_debug(logger, f'Replacing {(gaps > 0).sum()} gaps in the data with NaN entries')
    data = fill_gaps_cnsistent_wrapper(data)
    gaps_new = ((data.iloc[1:]['start'].values - data.iloc[:-1]['end'].values)
                * (data.iloc[:-1]['sample_id'].values == data.iloc[1:]['sample_id'].values))
    assert (gaps_new > 0).sum() == 0
    assert (np.sort(all_sample_chrom_ids) == data.groupby(['sample_id', 'chrom']).size().index.sort_values()).all()
    assert total_cn_sum == data[cn_columns].sum().sum()
    assert (data.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
    assert data.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'

    # Filter
    ## Filter based on coverage
    log_debug(logger, f"Filtering out chromosomes that have more than 75% NaN entries")
    coverage = 1 - ((data
        .set_index(['sample_id', 'chrom', 'start', 'end'])
        .isna().any(axis=1).to_frame('isna')
        .eval('(end-start)*isna')
        .groupby(['sample_id', 'chrom'])
        .sum()
        ) / 
        (data
        .set_index(['sample_id', 'chrom', 'start', 'end'])
        .isna().any(axis=1).to_frame('isna')
        .eval('(end-start)')
        .groupby(['sample_id', 'chrom'])
        .sum()
        )).fillna(0).to_frame('coverage')
    coverage = coverage.join(coverage.groupby('chrom').max().rename(columns={'coverage': 'max_coverage'}), on='chrom')
    log_debug(logger, f"Found {coverage.eval('coverage/max_coverage < 0.75').sum()} ({coverage.eval('coverage/max_coverage < 0.75').mean():.2f}%) chromosomes with less than 75% coverage")
    data = data.set_index(['sample_id', 'chrom']).loc[coverage.query('coverage/max_coverage >= 0.75').index].reset_index()

    ## Filter out diploid chromosomes
    # remove diploid and NaN chroms
    neutral_value = pd.DataFrame(2 if total_cn else 1, index=np.concatenate([wgd_samples, nowgd_samples]), columns=['neutral_value'])
    neutral_value.loc[wgd_samples] = 4 if total_cn else 2
    if total_cn:
        non_diploid_chroms = (data
        .set_index(['sample_id', 'chrom'])
        .join(neutral_value, on='sample_id')
        .query('not total_cn.isna()')
        .eval('total_cn != neutral_value', engine='python')
        .groupby(['sample_id', 'chrom'])
        .any()
        )
    else:
        non_diploid_chroms = (data
        .set_index(['sample_id', 'chrom'])
        .join(neutral_value, on='sample_id')
        .query('not cn_a.isna()')
        .eval('cn_a != neutral_value or cn_b != neutral_value', engine='python')
        .groupby(['sample_id', 'chrom'])
        .any()
        )
    data = data.join(non_diploid_chroms.to_frame('non_diploid'), on=['sample_id', 'chrom']).query('non_diploid').drop(columns='non_diploid')
    data = data.reset_index(drop=True)
    log_debug(logger, f"Selected {non_diploid_chroms.sum()} chromosomes out of {len(non_diploid_chroms)} ({100* non_diploid_chroms.sum() / len(non_diploid_chroms):.1f}%) that are non-diploid")
    all_sample_chrom_ids = data.groupby(['sample_id', 'chrom']).size().index.sort_values()

    # Phase data using MEDICC2 phasing
    if not skip_phasing and not total_cn:
        total_cn_sum = data[cn_columns].sum().sum()
        log_debug(logger, 'Phasing data')
        non_nan_index = ~data[cn_columns[0]].isna()
        data_ = data.copy()
        data_ = data_.loc[non_nan_index].set_index(['sample_id', 'chrom', 'start', 'end'])[['cn_a', 'cn_b']].astype(int).astype(str).copy()

        phased_data = []
        for idx, df in data_.groupby('sample_id'):

            cur_fst = T_forced_WGD if idx in wgd_samples else nowgd_fst
            phasing_dict = medicc.create_phasing_fsa_dict_from_df(df, cur_fst.input_symbols(), 'X')
            fsa_dict_a, fsa_dict_b, scores = medicc.phase_dict(phasing_dict, cur_fst, diploid_fsa)
            cur_phased_data = medicc.core.create_df_from_phasing_fsa(df, [fsa_dict_a, fsa_dict_b], 'X')
            cur_phased_data.columns = ['cn_a', 'cn_b']
            cur_phased_data[['cn_a', 'cn_b']] = cur_phased_data[['cn_a', 'cn_b']].astype(int)

            phasing_improves_chrom = []
            for chrom in df.index.get_level_values('chrom').unique():
                score = float(fstlib.score(cur_fst, diploid_fsa, create_standard_fsa_dict_from_data(df.query('chrom == @chrom'), SYMBOL_TABLE)[0][idx]))
                score_phased = float(fstlib.score(cur_fst, diploid_fsa, create_standard_fsa_dict_from_data(cur_phased_data.query('chrom == @chrom'), SYMBOL_TABLE)[0][idx]))
                if score_phased < score:
                    phasing_improves_chrom.append(chrom)
            if len(phasing_improves_chrom) == 0:
                phased_data.append(df)
                continue

            cur_phased_data = cur_phased_data.query('chrom in @phasing_improves_chrom')

            # flip groups of segments to always have cn_a "on top" if possible
            # neutral position is either cn_a == cn_b != 0 or chromosome break
            neutral_position = ((cur_phased_data['cn_a'] == cur_phased_data['cn_b']) & (cur_phased_data['cn_a'] != 0))
            neutral_position = neutral_position | np.append(False, (cur_phased_data.reset_index()['chrom'].values[:-1] != cur_phased_data.reset_index()['chrom'].values[1:]))

            group_boundaries = np.unique(np.concatenate([np.where(neutral_position)[0], [-1, len(cur_phased_data)]]))
            for i, j in zip(group_boundaries[:-1], group_boundaries[1:]):
                cur_vals = cur_phased_data.iloc[i:j][['cn_a', 'cn_b']].values.copy()
                if (cur_vals[:, 0] < cur_vals[:, 1]).sum() > (cur_vals[:, 0] > cur_vals[:, 1]).sum():
                    cur_phased_data.iloc[i:j] = cur_vals[:, [1, 0]]
            # Flip such that cn_b is zero for XY samples on chrX and chrY
            if idx in xy_samples:
                for cur_chrom in ['chrX', 'chrY']:
                    cur_ind = cur_phased_data.query('chrom == @cur_chrom').index
                    if ((cur_phased_data.query('chrom == @cur_chrom')['cn_a'] == 0).all() and
                        (cur_phased_data.query('chrom == @cur_chrom')['cn_b'] != 0).any()):
                        cur_phased_data.loc[cur_ind, 'cn_a'] = cur_phased_data.loc[cur_ind, 'cn_b']
                        cur_phased_data.loc[cur_ind, 'cn_b'] = 0

            cur_phased_data = cur_phased_data.astype(int)
            
            cur_phased_data = pd.concat(
                [df.query('chrom not in @phasing_improves_chrom'),
                cur_phased_data.query('chrom in @phasing_improves_chrom')]).sort_index()
            assert (df.astype(int).sum(axis=1) == cur_phased_data.astype(int).sum(axis=1)).all()
            
            if idx in xy_samples:
                cur_phased_data.loc[cur_phased_data.query('(chrom == "chrX") or (chrom == "chrY")').index, 'cn_a'] = df.loc[df.query('(chrom == "chrX") or (chrom == "chrY")').index, 'cn_a'].astype(int)
                cur_phased_data.loc[cur_phased_data.query('(chrom == "chrX") or (chrom == "chrY")').index, 'cn_b'] = 0

            phased_data.append(cur_phased_data)
        del(data_)
        phased_data = pd.concat(phased_data)
        data.loc[non_nan_index, ['cn_a', 'cn_b']] = phased_data[['cn_a', 'cn_b']].astype(int).values

        assert total_cn_sum == data[['cn_a', 'cn_b']].sum().sum()
        assert (all_sample_chrom_ids == data.groupby(['sample_id', 'chrom']).size().index.sort_values()).all()
        assert (data.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
        assert data.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'
    else:
        log_debug(logger, 'Skipping phasing')

    ## Throw out duplicate chromosomes
    if unique_chroms:
        log_debug(logger, 'Looking for duplicate chromosomes')

        alleles = ['cn_a'] if total_cn else ['cn_a', 'cn_b']
        duplicate_chrom_mapping = []
        for allele in alleles:
            for cur_wgd, cur_wgd_samples in zip(['noWGD', 'WGD'], [nowgd_samples, wgd_samples]):
                cn_strings = data.query('sample_id in @cur_wgd_samples').groupby(['sample_id', 'chrom'])[allele].agg(lambda x: ''.join(x.astype(str)))
                # add chrom here to make sure that it always matches the same chrom
                cn_strings.loc[:] = cn_strings.values + cn_strings.index.get_level_values('chrom').astype(str)
                log_debug(logger, f'For {cur_wgd} and allele {allele} found {len(cn_strings) - len(cn_strings.drop_duplicates())} ({(len(cn_strings) - len(cn_strings.drop_duplicates())) / len(cn_strings) * 100:.2f}%) duplicate IDs')

                cn_strings.index = cn_strings.index.map(lambda x: f'{x[0]}:{x[1]}')
                cn_strings_unique = cn_strings.drop_duplicates(keep='first')
                cn_strings_unique_dict = {v: k for k, v in cn_strings_unique.to_dict().items()}
                cur_duplicate_chrom_mapping = cn_strings.map(cn_strings_unique_dict)
                cur_duplicate_chrom_mapping = cur_duplicate_chrom_mapping.iloc[cur_duplicate_chrom_mapping.index.values != cur_duplicate_chrom_mapping.values]

                cur_duplicate_chrom_mapping.loc[:] = cur_duplicate_chrom_mapping.values + ':' + allele
                cur_duplicate_chrom_mapping.index = cur_duplicate_chrom_mapping.index.values + ':' + allele

                duplicate_chrom_mapping.append(cur_duplicate_chrom_mapping)
        duplicate_chrom_mapping = pd.concat(duplicate_chrom_mapping).reset_index().rename(columns={'index': 'id', 0: 'target'})
        duplicate_chrom_mapping.to_csv(os.path.join(directories['data_dir'], f'{name}_duplicate_chrom_mapping.tsv'), sep='\t', index=False)

        duplicate_chroms = pd.concat([
            pd.DataFrame(duplicate_chrom_mapping['id'].str.split(':', expand=True).rename(columns={0: 'sample_id', 1: 'chrom', 2: 'allele'})),
            pd.DataFrame(duplicate_chrom_mapping['target'].str.split(':', expand=True).rename(columns={0: 'target_sample_id', 1: 'target_chrom', 2: 'target_allele'}))

        ], axis=1)
        assert duplicate_chroms.eval('chrom == target_chrom').all(), f'duplicate chrom mapping maps between different chromosomes: {duplicate_chroms.query("chrom != target_chrom")}'
        assert duplicate_chrom_mapping.eval('id == target').sum() == 0, 'there are still some IDs left that map to themselves'
        log_debug(logger, f"Found a total of {duplicate_chrom_mapping['id'].nunique()} duplicate IDs that match to {duplicate_chrom_mapping['target'].nunique()} unique IDs")

    # Merge neighboring segments without cn change
    old_len = len(data)
    data = merge_neighbours_mod(
        data, cn_columns=cn_columns_phased, start_end_must_overlap=False)
    assert (data.reset_index().index == data.index).all(), 'index is messed up'
    assert (all_sample_chrom_ids == data.groupby(['sample_id', 'chrom']).size().index.sort_values()).all()
    assert (data.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
    assert data.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'

    log_debug(logger, f'Merged {old_len - len(data)} segments without copynumber change')

    # Telomeres, Centromeres and short arms
    ## Telomeres (up to 1Mb of next one)
    old_len = len(data)
    old_nan = data[cn_columns_phased[0]].isna().sum()
    data = fill_telomere_nans(data, cn_columns=cn_columns_phased)
    data.loc[data.query('(end - start) < 1e3').index.values, 'cn_a'] = np.nan
    log_debug(logger, f'Filled {old_nan - data[cn_columns_phased[0]].isna().sum()} telomeric NaN entries')

    assert data.query('(end - start) < 1e3')['cn_a'].isna().all()
    assert (all_sample_chrom_ids == data.groupby(['sample_id', 'chrom']).size().index.sort_values()).all()
    assert (data.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
    assert data.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'

    ## Unify centromeres
    if skip_centromeres:
        log_debug(logger, 'Skipping centromere binning')
    else:
        log_debug(logger, 'Creating uniform centromere bins and merging with data (might take a while)')
        centromeres = data_loaders.load_centromeres(extended=True)
        # assert that distances of centromeres to telomeres are bigger than 1e6 so the telomere NaN merging should be fine
        assert centromeres.join(pd.Series(chrom_lens_df, name='chrom_length'), on='chrom').astype(int).eval('chrom_length - centro_end').min() > 1e6

        centromeres_df = centromeres.reset_index().rename({'centro_start': 'start', 'centro_end': 'end'}, axis=1).copy()
        centromeres_df['sample_id'] = 'centro-' + centromeres_df['chrom']
        centromeres_df['cn'] = -1
        centromeres_df = centromeres_df[['sample_id', 'chrom', 'start', 'end', 'cn']]
        centromeres_df = fill_gaps_cnsistent_wrapper(centromeres_df, print_info=False)
        old_len = len(data)
        old_nan = data['cn_a'].isna().sum()
        for col, col_phased in zip(cn_columns, cn_columns_phased):
            data = data.rename(columns={col_phased: col})
        data['sample_chrom_id'] = data['sample_id'] + ':' + data['chrom']
        data[cn_columns] = data[cn_columns].astype(float)

        all_data = []
        for cur_id in data['sample_chrom_id'].unique():
            cur_chrom = cur_id.split(':')[1]
            cur_data = data.query('sample_chrom_id == @cur_id').reset_index(drop=True).copy().drop(columns='sample_chrom_id')
            if not cur_data[cn_columns[0]].isna().any():
                all_data.append(cur_data)
                continue

            breaks = get_breaks_mod(pd.concat([cur_data, centromeres_df.query('chrom == @cur_chrom')]))
            segments = cnsistent_breaks_to_segments(breaks)
            cur_data_with_centro = main_aggregate_quiet(cur_data, segments, cn_columns=cn_columns, how="mean", print_info=False)
            cur_data_with_centro = cur_data_with_centro.drop(columns='name')
            centromere_index = cur_data_with_centro.loc[cur_data_with_centro['start'] == centromeres.loc[cur_chrom, 'centro_start']].index[0]
            # check that the original segment was not inside the centromere
            assert cur_data_with_centro.loc[centromere_index, 'end'] - centromeres.loc[cur_chrom, 'centro_end'] < 2
            if centromere_index > 1 and (
                    cur_data_with_centro.at[centromere_index-1, 'end'] - cur_data_with_centro.at[centromere_index-1, 'start'] < 1e6):
                for cn_col in cn_columns:
                    cur_data_with_centro.loc[centromere_index-1, cn_col] = cur_data_with_centro.loc[centromere_index-2, cn_col]
            if centromere_index < len(cur_data) - 1 and (
                    cur_data_with_centro.at[centromere_index+1, 'end'] - cur_data_with_centro.at[centromere_index+1, 'start'] < 1e6):
                for cn_col in cn_columns:
                    cur_data_with_centro.loc[centromere_index+1, cn_col] = cur_data_with_centro.loc[centromere_index+2, cn_col]

            assert (cur_data_with_centro.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
            assert cur_data_with_centro.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'

            all_data.append(cur_data_with_centro)

        assert len(data['sample_chrom_id'].unique()) == len(all_data)
        data = pd.concat(all_data, axis=0).reset_index(drop=True)
        del(all_data)

        cur = data.join(centromeres, on='chrom')
        centromere_index = ((cur['start'] - cur['centro_start']).abs() < 10) & ((cur['end'] - cur['centro_end']).abs() < 10) & cur[cn_columns[0]].isna()
        del(cur)
        for cn_col in cn_columns:
            data.loc[centromere_index, cn_col] = -1

        data = merge_neighbours_mod(data, cn_columns=cn_columns, start_end_must_overlap=False)

        for col, col_phased in zip(cn_columns, cn_columns_phased):
            data = data.rename(columns={col: col_phased})
        for col in cn_columns_phased:
            data.loc[data.query('(end - start) < 1e3').index.values, col] = np.nan

        log_debug(logger, f'Filled {old_len - len(data)} entries with {old_nan - data[cn_columns_phased[0]].isna().sum()} NaNs')

        if not total_cn:
            assert ((data['cn_a'] ==-1) == (data['cn_b'] ==-1)).all()
            assert ((data['cn_a'].isna()) == (data['cn_b'].isna())).all()
        assert data.query('(end - start) < 1e3')[cn_columns_phased[0]].isna().all()
        assert (all_sample_chrom_ids == data.groupby(['sample_id', 'chrom']).size().index.sort_values()).all()
        assert (data.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
        assert data.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'

    # Short arms (13, 14, 15, 21, 22)
    short_arm_chroms = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']
    short_arm_nans = data.query('cn_a.isna() and (start == 0) and chrom in @short_arm_chroms', engine='python').index.values
    data.loc[short_arm_nans, 'cn_a'] = -1
    if not total_cn:
        data.loc[short_arm_nans, 'cn_b'] = -1
    log_debug(logger, f'Filled {len(short_arm_nans)} NaNs in short arms')

    assert (all_sample_chrom_ids == data.groupby(['sample_id', 'chrom']).size().index.sort_values()).all()
    assert (data.groupby(['sample_id', 'chrom'])['start'].min() == 0).all(), 'not all chroms start at zero'
    assert data.groupby(['sample_id', 'chrom'])['end'].max().to_frame('max_end').join(chrom_lens_df, on='chrom').eval('max_end == chrom_length').all(), 'not all chroms end at chrom length'


    # Remove all NaNs
    log_debug(logger, f"Removing {(data['cn_a'] == -1).sum()} NaNs already flagged as centromeres and telomeres and {data['cn_a'].isna().sum()} NaNs not flagged yet")
    data = data.loc[(data['cn_a']!= -1) & ~(data['cn_a'].isna())]
    assert not data.isna().any().any(), data.isna().sum()
    assert not (data[cn_columns_phased] == -1).any().any()
    data = data.reset_index(drop=True)

    ## Final merge after NaN removal
    data = merge_neighbours_mod(data, cn_columns=cn_columns_phased, start_end_must_overlap=False)
    data = data.reset_index(drop=True)
    assert (np.sort(all_sample_chrom_ids) == np.sort(data.groupby(['sample_id', 'chrom']).size().index)).all()
    assert data.eval('end-start').min() > 1e3

    # Separate the two haplotypes
    if total_cn:
        data_stacked = data.copy()
        data_stacked = data_stacked.rename(columns={'cn_a': 'cn'})
        data_stacked['allele'] = 'cn_a'
    else:
        data.columns.name = 'allele'
        data_stacked = data.set_index(['sample_id', 'chrom', 'start', 'end']).stack().to_frame('cn').reset_index().copy()
        assert len(data) == len(data_stacked)//2, (len(data), len(data_stacked), len(data_stacked)//2)
        assert data_stacked['cn'].sum() == data[cn_columns_phased].sum().sum()
    data_stacked['cn'] = data_stacked['cn'].astype(int)
    data_stacked['id'] = data_stacked['sample_id'] + ':' + data_stacked['chrom'] + ':' + data_stacked['allele']

    ## Remove diploid chroms
    neutral_value = pd.DataFrame(2 if total_cn else 1, index=np.concatenate([wgd_samples, nowgd_samples]), columns=['neutral_value'])
    neutral_value.loc[wgd_samples] = 4 if total_cn else 2
    old_n_ids = data_stacked['id'].nunique()
    non_diploid_ids = (data_stacked
     .set_index(['id'])
     .join(neutral_value, on='sample_id')
     .eval('cn != neutral_value')
     .groupby(['id'])
     .any()
     .to_frame('non_diploid').query('non_diploid').index
     )
    data_stacked = data_stacked.query('id in @non_diploid_ids').reset_index(drop=True).copy()
    assert (data['start'] != 1).all()
    log_debug(logger, f"Selected {data_stacked['id'].nunique()} non-diploid IDs out of {old_n_ids} IDs")

    ## For total CN: Remove chrX and chrY with CN = 1 for male samples
    if total_cn:
        neutral_value = pd.DataFrame(1, index=np.concatenate([wgd_samples, nowgd_samples]), columns=['neutral_value'])
        neutral_value.loc[wgd_samples] = 2
        diploid_chrx_chry = (
            data_stacked
            .query('(chrom == "chrX" or chrom == "chrY") and sample_id in @xy_samples')
            .set_index(['id'])
            .join(neutral_value, on='sample_id')
            .eval('cn == neutral_value')
            .groupby(['id'])
            .all()
            .to_frame('diploid').query('diploid').index.values
        )
        data_stacked = data_stacked.query('id not in @diploid_chrx_chry').reset_index(drop=True).copy()
        log_debug(logger, f"Removed {len(diploid_chrx_chry)} chrX or chrY chroms for male samples with CN = 1")
    ## Remove Y chroms for female (XX) samples
    if xy_samples is not None:
        old_n_ids = data_stacked['id'].nunique()
        data_stacked = data_stacked.query('not ((sample_id not in @xy_samples) and chrom == "chrY")').copy().reset_index(drop=True)
        log_debug(logger, f"Removed {old_n_ids - data_stacked['id'].nunique()} chrY chroms from female (XX) samples")
        ## Remove minor X and Y chroms from male (XY) samples
        old_n_ids = data_stacked['id'].nunique()
        query = '(sample_id in @xy_samples) and allele == "cn_b" and (chrom == "chrX" or chrom == "chrY")'
        assert (data_stacked.query(query)['cn'] == 0).all()
        data_stacked = data_stacked.query(f'not ({query})').copy().reset_index(drop=True)
        ## For total CN: remove chrY from female (XX) samples
        if total_cn:
            old_n_ids = data_stacked['id'].nunique()
            data_stacked = data_stacked.query('not ((sample_id not in @xy_samples) and chrom == "chrY")').copy().reset_index(drop=True)
            log_debug(logger, f"Removed {old_n_ids - data_stacked['id'].nunique()} chrY chroms from female (XX) samples")

    ## Now remove all duplicate chroms
    if unique_chroms:
        duplicate_ids_to_remove = duplicate_chrom_mapping['id'].values
        old_n_ids = data_stacked['id'].nunique()
        data_stacked = data_stacked.query('id not in @duplicate_ids_to_remove').copy().reset_index(drop=True)
        new_n_ids = data_stacked['id'].nunique()
        log_debug(logger, f"Out of {old_n_ids} IDs, {old_n_ids - new_n_ids} duplicate chroms were removed leading to a final number of {new_n_ids} IDs")

    ##  Final merge (this is required for downstream analysis such as LOH adjust in full paths)
    data_stacked = data_stacked.sort_values(['id', 'start']).reset_index(drop=True)
    data_stacked['sample_id'] = data_stacked['sample_id'] + '_' + data_stacked['allele']
    data_stacked = merge_neighbours_mod(
        data_stacked, cn_columns=['cn'], start_end_must_overlap=False)
    data_stacked['sample_id'] = data_stacked['sample_id'].str.replace('_cn_a', '').str.replace('_cn_b', '')

    # Save final
    log_debug(logger, f'Save final to dir {directories["data_dir"]}')
    log_debug(logger, f'Final data: {name}_processed.tsv and {name}_processed_split.tsv')
    data_stacked = data_stacked.sort_values(['sample_id', 'chrom', 'allele'])
    data_stacked['id'] = data_stacked['sample_id'] + ':' + data_stacked['chrom'] + ':' + data_stacked['allele']
    data_stacked.to_csv(os.path.join(directories['data_dir'], f'{name}_processed_split.tsv'), sep='\t', index=False)
    data.to_csv(os.path.join(directories['data_dir'], f'{name}_processed.tsv'), sep='\t', index=False)

    log_debug(logger, 'Done.')

    return data, data_stacked


def cnsistent_breaks_to_segments(breakpoints):
    '''Copied from CNSistent v0.9.0'''
    segs = {}
    for chrom, breaks in breakpoints.items():
        segs[chrom] = []
        for i in range(len(breaks) - 1):
            segs[chrom].append((breaks[i], breaks[i + 1], f"{chrom}_{i}"))
    return segs