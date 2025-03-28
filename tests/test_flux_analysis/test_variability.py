"""Test functionalities of Flux Variability Analysis."""

import os
from typing import Callable, List

import numpy as np
import pandas as pd
import pytest

from cobra import Model
from cobra.exceptions import Infeasible
from cobra.flux_analysis.variability import (
    find_blocked_reactions,
    find_essential_genes,
    find_essential_reactions,
    flux_variability_analysis,
)


# FVA
def test_flux_variability_benchmark(
    large_model: Model, benchmark: Callable, all_solvers: List[str]
) -> None:
    """Benchmark FVA."""
    large_model.solver = all_solvers
    benchmark(
        flux_variability_analysis,
        large_model,
        reaction_list=large_model.reactions[1::100],
        processes=1,
    )


def test_pfba_flux_variability(
    model: Model,
    pfba_fva_results: pd.DataFrame,
    fva_results: pd.DataFrame,
    all_solvers: List[str],
) -> None:
    """Test FVA using pFBA."""
    model.solver = all_solvers
    with pytest.warns(UserWarning):
        flux_variability_analysis(
            model,
            pfba_factor=0.1,
            reaction_list=model.reactions[1::3],
            processes=1,
        )
    fva_out = flux_variability_analysis(
        model, pfba_factor=1.1, reaction_list=model.reactions, processes=1
    )
    fva_out.sort_index(inplace=True)
    assert np.allclose(fva_out, pfba_fva_results)
    abs_fva_out = fva_out.dropna().abs()
    abs_fva_results = fva_results.dropna().abs()
    comparison = np.isclose(abs_fva_out, abs_fva_results) | (
        abs_fva_out < abs_fva_results
    )
    assert comparison["minimum"].all()
    assert comparison["maximum"].all()


def test_loopless_pfba_fva(model: Model) -> None:
    """Test loopless FVA using pFBA."""
    loop_reactions = [model.reactions.get_by_id(rid) for rid in ("FRD7", "SUCDi")]
    fva_loopless = flux_variability_analysis(
        model,
        pfba_factor=1.1,
        reaction_list=loop_reactions,
        loopless=True,
        processes=1,
    )
    assert np.allclose(fva_loopless["maximum"], fva_loopless["minimum"])


def test_flux_variability(
    model: Model, fva_results: pd.DataFrame, all_solvers: List[str]
) -> None:
    """Test FVA."""
    model.solver = all_solvers
    fva_out = flux_variability_analysis(
        model, reaction_list=model.reactions, processes=1
    )
    fva_out.sort_index(inplace=True)
    assert np.allclose(fva_out, fva_results)


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_parallel_flux_variability(
    model: Model, fva_results: pd.DataFrame, all_solvers: List[str]
) -> None:
    """Test parallel FVA."""
    model.solver = all_solvers
    fva_out = flux_variability_analysis(model, processes=2)
    fva_out.sort_index(inplace=True)
    assert np.allclose(fva_out, fva_results)


# Loopless FVA
def test_flux_variability_loopless_benchmark(
    model: Model, benchmark: Callable, all_solvers: List[str]
) -> None:
    """Benchmark loopless FVA."""
    model.solver = all_solvers
    benchmark(
        flux_variability_analysis,
        model,
        loopless=True,
        reaction_list=model.reactions[1::50],
    )


def test_flux_variability_loopless(model: Model, all_solvers: List[str]) -> None:
    """Test loopless FVA."""
    model.solver = all_solvers
    loop_reactions = [model.reactions.get_by_id(rid) for rid in ("FRD7", "SUCDi")]
    fva_normal = flux_variability_analysis(model, reaction_list=loop_reactions)
    fva_loopless = flux_variability_analysis(
        model, reaction_list=loop_reactions, loopless=True
    )

    assert not np.allclose(fva_normal["maximum"], fva_normal["minimum"])
    assert np.allclose(fva_loopless["maximum"], fva_loopless["minimum"])


# Internals (essentiality, blocked reactions)
def test_fva_data_frame(model: Model) -> None:
    """Test DataFrame obtained from FVA."""
    df = flux_variability_analysis(model)
    assert np.all([df.columns.values == ["minimum", "maximum"]])


def test_fva_infeasible(model: Model) -> None:
    """Test FVA infeasibility."""
    infeasible_model = model.copy()
    infeasible_model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
    # ensure that an infeasible model does not run FVA
    with pytest.raises(Infeasible):
        flux_variability_analysis(infeasible_model)


def test_fva_minimization(model: Model) -> None:
    """Test minimization using FVA."""
    model.objective = model.reactions.EX_glc__D_e
    model.objective_direction = "min"
    solution = flux_variability_analysis(model, fraction_of_optimum=0.95)
    assert solution.at["EX_glc__D_e", "minimum"] == -10.0
    assert solution.at["EX_glc__D_e", "maximum"] == -9.5


def test_find_blocked_reactions_solver_none(model: Model) -> None:
    """Test find_blocked_reactions() [no specific solver]."""
    result = find_blocked_reactions(model, model.reactions[40:46])
    assert result == ["FRUpts2"]


def test_essential_genes(model: Model) -> None:
    """Test find_essential_genes()."""
    essential_genes = {
        "b2779",
        "b1779",
        "b0720",
        "b2416",
        "b2926",
        "b1136",
        "b2415",
    }
    observed_essential_genes = {g.id for g in find_essential_genes(model)}
    assert observed_essential_genes == essential_genes


def test_essential_reactions(model: Model) -> None:
    """Test find_blocked_reactions()."""
    essential_reactions = {
        "GLNS",
        "Biomass_Ecoli_core",
        "PIt2r",
        "GAPD",
        "ACONTb",
        "EX_nh4_e",
        "ENO",
        "EX_h_e",
        "EX_glc__D_e",
        "ICDHyr",
        "CS",
        "NH4t",
        "GLCpts",
        "PGM",
        "EX_pi_e",
        "PGK",
        "RPI",
        "ACONTa",
    }
    observed_essential_reactions = {r.id for r in find_essential_reactions(model)}
    assert observed_essential_reactions == essential_reactions


def test_find_blocked_reactions(model: Model, all_solvers: List[str]) -> None:
    """Test find_blocked_reactions()."""
    model.solver = all_solvers
    result = find_blocked_reactions(model, model.reactions[40:46])
    assert result == ["FRUpts2"]

    result = find_blocked_reactions(model, model.reactions[42:48])
    assert set(result) == {"FUMt2_2", "FRUpts2"}

    result = find_blocked_reactions(model, model.reactions[30:50], open_exchanges=True)
    assert result == []
