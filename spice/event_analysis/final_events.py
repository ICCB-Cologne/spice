def classify_event_position(final_events_df):
    final_events_df['pos'] = final_events_df[['whole_chrom', 'whole_arm', 'telomere_bound']].apply(
        lambda x: 'whole_chrom' if x['whole_chrom'] else (
            'whole_arm' if x['whole_arm'] else (
                'telomere_bound' if x['telomere_bound'] else 'internal'
            )
        ), axis=1
    )
    return final_events_df['pos']
