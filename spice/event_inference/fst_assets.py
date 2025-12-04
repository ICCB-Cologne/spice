import medicc
from spice.event_inference.fsts import *


max_pre_wgd_losses = 8
separator = 'X'
max_cn = 8
SYMBOL_TABLE = medicc.create_symbol_table(max_cn=max_cn, separator=separator)

# diploid_fsa = medicc.tools.create_diploid_fsa(medicc.io.read_fst())
def get_diploid_fsa(total_copy_numbers=False):
    def _create_diploid_fsa(fst, total_copy_numbers=False):
        '''Copied from medicc.tools.create_diploid_fsa'''
        import fstlib

        # Create diploid FSA
        diploid_fsa = fstlib.Fst()
        diploid_fsa.set_input_symbols(fst.input_symbols())
        diploid_fsa.set_output_symbols(fst.output_symbols())
        diploid_fsa.add_state()
        diploid_fsa.set_start(0)
        diploid_fsa.set_final(0, 0)
        if total_copy_numbers:
            diploid_fsa.add_arc(0, ('2', '2', 0, 0))
        else:
            diploid_fsa.add_arc(0, ('1', '1', 0, 0))
        diploid_fsa.add_arc(0, ('X', 'X', 0, 0))

        return diploid_fsa
    return _create_diploid_fsa(medicc.io.read_fst(), total_copy_numbers=total_copy_numbers)

# Create usefull FSTs
n = len(medicc.factory._get_int_cns_from_symbol_table(SYMBOL_TABLE, separator))
L1step = medicc.factory.create_1step_del_fst(SYMBOL_TABLE, separator, exclude_zero=True)
G1step = ~L1step
G = medicc.create_nstep_fst(n-1, G1step)
L = ~G
G_pre = medicc.create_nstep_fst(3, G1step)
LOH = medicc.factory.create_loh_fst(SYMBOL_TABLE)
L_LOH_1step = medicc.create_1step_del_fst(SYMBOL_TABLE, separator, exclude_zero=False)
L_LOH = medicc.create_nstep_fst(max_pre_wgd_losses-1, L_LOH_1step)
XX = fstlib.encode_determinize_minimize(L_LOH*G)
XX_noLOH = fstlib.encode_determinize_minimize(L*G)

W_forced = create_1step_forced_WGD_fst(
    SYMBOL_TABLE, separator, wgd_cost=1, minimize=False, wgd_x2=True, total_cn=False)
W = medicc.factory.create_1step_WGD_fst(
    SYMBOL_TABLE, separator, wgd_cost=1, minimize=False, wgd_x2=True, total_cn=False)

T_forced_WGD = L_LOH * G_pre * W_forced * XX
T_forced_WGD_without_post_LOH = LOH * G_pre * W_forced * XX_noLOH
T = LOH * G_pre * W * XX
T_noWGD = LOH * XX

wgd_fst = medicc.io.read_fst()
nowgd_fst = medicc.io.read_fst(no_wgd=True)
