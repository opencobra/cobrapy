""" Test functionalities of find_active_reactions"""

from __future__ import absolute_import

import pytest

import cobra
from cobra.test import create_test_model
from cobra.flux_analysis.find_active_reactions import (
    find_active_reactions, find_reactions_in_cycles)


def test_find_active_reactions_benchmark(model, benchmark, all_solvers):
    """Benchmark find_active_reactions."""

    model.solver = all_solvers
    benchmark(find_active_reactions, model)


def test_find_reactions_in_cycles_benchmark(model, benchmark, all_solvers):
    """Benchmark find_reactions_in_cycles."""

    model.solver = all_solvers
    benchmark(find_reactions_in_cycles, model)


# def test_fastSNP_benchmark(model, benchmark, all_solvers):
#     """Benchmark fastSNP."""
#
#     model.solver = all_solvers
#     benchmark(fastSNP, model)


def test_find_active_reactions(model, all_solvers):
    """Test find_active_reactions."""

    model.solver = all_solvers
    # solve LPs
    active_rxns_lp = find_active_reactions(model)

    active_rxns = ['ACALD', 'ACALDt', 'ACKr', 'ACONTa', 'ACONTb', 'ACt2r',
                   'ADK1', 'AKGDH', 'AKGt2r', 'ALCD2x', 'ATPM', 'ATPS4r',
                   'Biomass_Ecoli_core', 'CO2t', 'CS', 'CYTBD',
                   'D_LACt2', 'ENO', 'ETOHt2r', 'EX_ac_e', 'EX_acald_e',
                   'EX_akg_e', 'EX_co2_e', 'EX_etoh_e', 'EX_for_e',
                   'EX_glc__D_e', 'EX_glu__L_e', 'EX_h_e', 'EX_h2o_e',
                   'EX_lac__D_e', 'EX_nh4_e', 'EX_o2_e', 'EX_pi_e',
                   'EX_pyr_e', 'EX_succ_e', 'FBA', 'FBP', 'FORt2', 'FORti',
                   'FRD7', 'FUM', 'G6PDH2r', 'GAPD', 'GLCpts', 'GLNS',
                   'GLUDy', 'GLUN', 'GLUSy', 'GLUt2r', 'GND', 'H2Ot',
                   'ICDHyr', 'ICL', 'LDH_D', 'MALS', 'MDH', 'ME1', 'ME2',
                   'NADH16', 'NADTRHD', 'NH4t', 'O2t', 'PDH', 'PFK', 'PFL',
                   'PGI', 'PGK', 'PGL', 'PGM', 'PIt2r', 'PPC', 'PPCK', 'PPS',
                   'PTAr', 'PYK', 'PYRt2', 'RPE', 'RPI', 'SUCCt2_2', 'SUCCt3',
                   'SUCDi', 'SUCOAS', 'TALA', 'THD2', 'TKT1', 'TKT2', 'TPI']

    assert set(active_rxns_lp) == set(active_rxns)

    # solving MILP or Fast-SNP may not be feasible for some solvers
    if all_solvers in ["gurobi", "cplex"]:
        active_rxns_milp = find_active_reactions(model, solve="milp")
        active_rxns_fastsnp = find_active_reactions(model, solve="fastSNP")
        assert set(active_rxns_milp) == set(active_rxns)
        assert set(active_rxns_fastsnp) == set(active_rxns)


def test_find_reactions_in_cycles(large_model, all_solvers):
    """Test find_reactions_in_cycles."""

    large_model.solver = all_solvers
    # solve LPs
    rxns_in_cycles_lp = find_reactions_in_cycles(large_model)

    rxns_in_cycles = ['ABUTt2pp', 'ACCOAL', 'ACKr', 'ACS', 'ACt2rpp',
                      'ACt4pp', 'ADK1', 'ADK3', 'ADNt2pp', 'ADNt2rpp',
                      'ALATA_L', 'ALAt2pp', 'ALAt2rpp', 'ALAt4pp', 'ASPt2pp',
                      'ASPt2rpp', 'CA2t3pp', 'CAt6pp', 'CRNDt2rpp',
                      'CRNt2rpp', 'CRNt8pp', 'CYTDt2pp', 'CYTDt2rpp',
                      'FOMETRi', 'GLBRAN2', 'GLCP', 'GLCP2', 'GLCS1',
                      'GLCtex', 'GLCtexi', 'GLDBRAN2', 'GLGC', 'GLUABUTt7pp',
                      'GLUt2rpp', 'GLUt4pp', 'GLYCLTt2rpp', 'GLYCLTt4pp',
                      'GLYt2pp', 'GLYt2rpp', 'GLYt4pp', 'HPYRI', 'HPYRRx',
                      'ICHORS', 'ICHORSi', 'INDOLEt2pp', 'INDOLEt2rpp',
                      'INSt2pp', 'INSt2rpp', 'NAt3pp', 'NDPK1', 'PPAKr',
                      'PPCSCT', 'PPKr', 'PPM', 'PROt2rpp', 'PROt4pp',
                      'PRPPS', 'PTA2', 'PTAr', 'R15BPK', 'R1PK', 'SERt2rpp',
                      'SERt4pp', 'SUCOAS', 'THFAT', 'THMDt2pp', 'THMDt2rpp',
                      'THRt2rpp', 'THRt4pp', 'TRSARr', 'URAt2pp', 'URAt2rpp',
                      'URIt2pp', 'URIt2rpp', 'VALTA', 'VPAMTr']

    assert set(rxns_in_cycles_lp) == set(rxns_in_cycles)

    # solving MILP or Fast-SNP may not be feasible for some solvers
    if all_solvers in ["gurobi", "cplex"]:
        rxns_in_cycles_milp = find_reactions_in_cycles(large_model,
                                                       solve="milp")
        rxns_in_cycles_fastsnp = find_reactions_in_cycles(large_model,
                                                          solve="fastSNP")
        assert set(rxns_in_cycles_milp) == set(rxns_in_cycles)
        assert set(rxns_in_cycles_fastsnp) == set(rxns_in_cycles)
