import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from spice import data_loaders
from spice.event_analysis.final_events import classify_event_position

CHROM_LENS = data_loaders.load_chrom_lengths()
CHROMS = list(CHROM_LENS.keys())
CENTROMERES = data_loaders.load_centromeres()

def create_segmentation(size):
    cur_breakpoint_dict = {
        chrom: np.append((np.arange(0, CHROM_LENS[chrom], size))[:-1], CHROM_LENS[chrom])
        for chrom in CHROMS}

    cur_segmentation = pd.DataFrame(
        index=pd.MultiIndex.from_tuples(
            [(chrom, cur_breakpoint_dict[chrom][i], cur_breakpoint_dict[chrom][i+1]-1) 
             for chrom in CHROMS for i in range(len(cur_breakpoint_dict[chrom])-1)],
             names=['chrom', 'start', 'end']))
    
    return cur_segmentation

def create_events_in_segmentation(final_events_df, segmentation=100e3, show_tqdm=False):
    """
    Map events to segmentation bins.

    Returns a DataFrame indexed like `segmentation.index` with a single column
    `'events'` containing lists of event indices from `final_events_df` that
    overlap each segmentation bin.
    """
    if hasattr(segmentation, '__int__'):
        segmentation = create_segmentation(segmentation)
    else:
        assert isinstance(segmentation, pd.DataFrame)
    if 'pos' not in final_events_df.columns:
        final_events_df['pos'] = classify_event_position(final_events_df)

    events_in_segmentation = pd.DataFrame(
        0,
        index=segmentation.index,
        columns=['events']).sort_index()

    for cur_chrom in tqdm(final_events_df['chrom'].unique(), disable=not show_tqdm):
        cur_events = final_events_df.loc[final_events_df['chrom'] == cur_chrom]
        if len(cur_events) == 0:
            continue
        events_starts = cur_events['start'].values
        events_ends = cur_events['end'].values
        cur_seg = segmentation.loc[[cur_chrom]].reset_index().copy()
        cur_seg['end'] += 1

        cur_breakpoints = np.unique(cur_seg[['start', 'end']].values.flat)
        events_starts_smaller = events_starts < cur_breakpoints[1:, None]
        events_ends_larger = events_ends > cur_breakpoints[:-1, None]

        output = np.matmul(np.logical_and(events_starts_smaller, events_ends_larger), np.ones(len(cur_events), dtype=int))
        events_in_segmentation.loc[cur_chrom, 'events'] = output

    return events_in_segmentation

def create_events_in_segmentation_full(final_events_df, segmentation=100e3, show_tqdm=False):
    if hasattr(segmentation, '__int__'):
        segmentation = create_segmentation(segmentation)
    else:
        assert isinstance(segmentation, pd.DataFrame)
    if 'pos' not in final_events_df.columns:
        final_events_df['pos'] = classify_event_position(final_events_df)

    all_segmented_events = []
    for cur_type in ['gain', 'loss']:
        for cur_pos in ['whole_chrom', 'whole_arm', 'telomere_bound', 'internal']:
            cur_events = final_events_df.loc[(final_events_df['type'] == cur_type) & (final_events_df['pos'] == cur_pos)]
            cur_segmented_events = create_events_in_segmentation(cur_events, segmentation)
            all_segmented_events.append(cur_segmented_events)
    combined_segmented_events = pd.concat(all_segmented_events, axis=1)
    combined_segmented_events.columns = pd.MultiIndex.from_product(
        [['gain', 'loss'], ['whole_chrom', 'whole_arm', 'telomere_bound', 'internal']],
        names=['type', 'pos'])
    return combined_segmented_events


