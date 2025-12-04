from collections import namedtuple

ChromData = namedtuple(
    'ChromData',
    ['id', 'sample', 'chrom', 'allele', 'cn_profile', 'string', 'dist', 'n_events', 'has_wgd', 'copynumber_file'])
Diff = namedtuple(
    'Diff',
    ['diff', 'is_gain', 'wgd'])
FullPaths = namedtuple(
    'FullPaths', 
    ['id', 'sample', 'chrom', 'allele', 'cn_profile', 'n_solutions', 'n_events', 'is_wgd', 'solved', 'events',
     'solutions'])
