# -*- coding: utf-8 -*-

"""Test functionalities of reaction and gene deletions."""

from __future__ import absolute_import

import math

import numpy as np
import pytest
from pandas import Series
from six import iteritems

from cobra.flux_analysis.deletion import (
    double_gene_deletion,
    double_reaction_deletion,
    single_gene_deletion,
    single_reaction_deletion,
)
from cobra.flux_analysis.room import add_room


# Single gene deletion FBA
def test_single_gene_deletion_fba_benchmark(model, benchmark, all_solvers):
    """Benchmark single gene deletion using FBA."""
    model.solver = all_solvers
    benchmark(single_gene_deletion, model)


def test_single_gene_deletion_fba(model, all_solvers):
    """Test single gene deletion using FBA."""
    # expected knockouts for textbook model
    model.solver = all_solvers
    growth_dict = {
        "b0008": 0.87,
        "b0114": 0.80,
        "b0116": 0.78,
        "b2276": 0.21,
        "b1779": 0.00,
    }
    result = single_gene_deletion(
        model=model, gene_list=list(growth_dict), method="fba", processes=1
    )
    for gene, value in iteritems(growth_dict):
        assert np.isclose(result.knockout[gene].growth, value, atol=1e-02)


# Singe gene deletion MOMA
def test_single_gene_deletion_moma_benchmark(model, benchmark, qp_solvers):
    """Benchmark single gene deletion using MOMA."""
    model.solver = qp_solvers
    genes = ["b0008", "b0114", "b2276", "b1779"]
    benchmark(
        single_gene_deletion, model=model, gene_list=genes, method="moma", processes=1
    )


def test_single_gene_deletion_moma(model, qp_solvers):
    """Test single gene deletion using MOMA."""
    model.solver = qp_solvers
    # expected knockouts for textbook model
    growth_dict = {
        "b0008": 0.87,
        "b0114": 0.71,
        "b0116": 0.56,
        "b2276": 0.11,
        "b1779": 0.00,
    }

    result = single_gene_deletion(
        model=model, gene_list=list(growth_dict), method="moma", processes=1
    )
    for gene, value in iteritems(growth_dict):
        assert np.isclose(result.knockout[gene].growth, value, atol=1e-02)


def test_single_gene_deletion_moma_reference(model, qp_solvers):
    """Test single gene deletion using MOMA (reference solution)."""
    model.solver = qp_solvers
    # expected knockouts for textbook model
    growth_dict = {
        "b0008": 0.87,
        "b0114": 0.71,
        "b0116": 0.56,
        "b2276": 0.11,
        "b1779": 0.00,
    }

    sol = model.optimize()
    result = single_gene_deletion(
        model=model,
        gene_list=list(growth_dict),
        method="moma",
        solution=sol,
        processes=1,
    )
    for gene, value in iteritems(growth_dict):
        assert np.isclose(result.knockout[gene].growth, value, atol=1e-02)


# Single gene deletion linear MOMA
def test_single_gene_deletion_linear_moma_benchmark(model, benchmark, all_solvers):
    """Benchmark single gene deletion using linear MOMA."""
    model.solver = all_solvers
    genes = ["b0008", "b0114", "b2276", "b1779"]
    benchmark(
        single_gene_deletion,
        model=model,
        gene_list=genes,
        method="linear moma",
        processes=1,
    )


def test_single_gene_deletion_linear_moma(model, all_solvers):
    """Test single gene deletion using linear MOMA (reference solution)."""
    model.solver = all_solvers
    # expected knockouts for textbook model
    growth_dict = {
        "b0008": 0.87,
        "b0114": 0.76,
        "b0116": 0.65,
        "b2276": 0.08,
        "b1779": 0.00,
    }

    sol = model.optimize()
    result = single_gene_deletion(
        model=model,
        gene_list=list(growth_dict),
        method="linear moma",
        solution=sol,
        processes=1,
    )
    for gene, value in iteritems(growth_dict):
        assert np.isclose(result.knockout[gene].growth, value, atol=1e-02)


# Single gene deletion ROOM
def test_single_gene_deletion_room_benchmark(model, benchmark, all_solvers):
    """Benchmark single gene deletion using ROOM."""
    if all_solvers == "glpk":
        pytest.skip("GLPK is too slow to run ROOM.")
    model.solver = all_solvers
    genes = ["b0008", "b0114", "b2276", "b1779"]
    benchmark(
        single_gene_deletion, model=model, gene_list=genes, method="room", processes=1
    )


# Single gene deletion linear ROOM
def test_single_gene_deletion_linear_room_benchmark(model, benchmark, all_solvers):
    """Benchmark single gene deletion using linear ROOM."""
    model.solver = all_solvers
    genes = ["b0008", "b0114", "b2276", "b1779"]
    benchmark(
        single_gene_deletion,
        model=model,
        gene_list=genes,
        method="linear room",
        processes=1,
    )


# Single reaction deletion
def test_single_reaction_deletion_benchmark(model, benchmark, all_solvers):
    """Benchmark single reaction deletion."""
    model.solver = all_solvers
    benchmark(single_reaction_deletion, model=model, processes=1)


def test_single_reaction_deletion(model, all_solvers):
    """Test single reaction deletion."""
    model.solver = all_solvers
    expected_results = {
        "FBA": 0.70404,
        "FBP": 0.87392,
        "CS": 0,
        "FUM": 0.81430,
        "GAPD": 0,
        "GLUDy": 0.85139,
    }
    result = single_reaction_deletion(
        model=model, reaction_list=list(expected_results), processes=1
    )

    for reaction, value in iteritems(expected_results):
        assert np.isclose(result.knockout[reaction].growth, value, atol=1e-05)


# Single reaction deletion ROOM
def test_single_reaction_deletion_room(room_model, room_solution, all_solvers):
    """Test single reaction deletion using ROOM."""
    room_model.solver = all_solvers
    expected = Series(
        {
            "v1": 10.0,
            "v2": 5.0,
            "v3": 0.0,
            "v4": 5.0,
            "v5": 5.0,
            "v6": 0.0,
            "b1": 10.0,
            "b2": 5.0,
            "b3": 5.0,
        },
        index=["v1", "v2", "v3", "v4", "v5", "v6", "b1", "b2", "b3"],
    )
    with room_model:
        room_model.reactions.v6.knock_out()
        add_room(room_model, solution=room_solution, delta=0.0, epsilon=0.0)
        room_sol = room_model.optimize()

    assert np.allclose(room_sol.fluxes, expected)


# Single reaction deletion linear ROOM
def test_single_reaction_deletion_linear_room(room_model, room_solution, all_solvers):
    """Test single reaction deletion using linear ROOM."""
    room_model.solver = all_solvers
    expected = Series(
        {
            "v1": 10.0,
            "v2": 5.0,
            "v3": 0.0,
            "v4": 5.0,
            "v5": 5.0,
            "v6": 0.0,
            "b1": 10.0,
            "b2": 5.0,
            "b3": 5.0,
        },
        index=["v1", "v2", "v3", "v4", "v5", "v6", "b1", "b2", "b3"],
    )
    with room_model:
        room_model.reactions.v6.knock_out()
        add_room(
            room_model, solution=room_solution, delta=0.0, epsilon=0.0, linear=True
        )
        linear_room_sol = room_model.optimize()

    assert np.allclose(linear_room_sol.fluxes, expected)


# Double gene deletion
def test_double_gene_deletion_benchmark(large_model, benchmark):
    """Benchmark double gene deletion."""
    genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935", "b1276", "b1241"]
    benchmark(double_gene_deletion, large_model, gene_list1=genes, processes=1)


def test_double_gene_deletion(model):
    """Test double gene deletion."""
    genes = ["b0726", "b4025", "b0724", "b0720", "b2935", "b2935", "b1276", "b1241"]
    growth_dict = {
        "b0720": {
            "b0720": 0.0,
            "b0724": 0.0,
            "b0726": 0.0,
            "b1241": 0.0,
            "b1276": 0.0,
            "b2935": 0.0,
            "b4025": 0.0,
        },
        "b0724": {
            "b0720": 0.0,
            "b0724": 0.814,
            "b0726": 0.814,
            "b1241": 0.814,
            "b1276": 0.814,
            "b2935": 0.814,
            "b4025": 0.739,
        },
        "b0726": {
            "b0720": 0.0,
            "b0724": 0.814,
            "b0726": 0.858,
            "b1241": 0.858,
            "b1276": 0.858,
            "b2935": 0.858,
            "b4025": 0.857,
        },
        "b1241": {
            "b0720": 0.0,
            "b0724": 0.814,
            "b0726": 0.858,
            "b1241": 0.874,
            "b1276": 0.874,
            "b2935": 0.874,
            "b4025": 0.863,
        },
        "b1276": {
            "b0720": 0.0,
            "b0724": 0.814,
            "b0726": 0.858,
            "b1241": 0.874,
            "b1276": 0.874,
            "b2935": 0.874,
            "b4025": 0.863,
        },
        "b2935": {
            "b0720": 0.0,
            "b0724": 0.814,
            "b0726": 0.858,
            "b1241": 0.874,
            "b1276": 0.874,
            "b2935": 0.874,
            "b4025": 0.863,
        },
        "b4025": {
            "b0720": 0.0,
            "b0724": 0.739,
            "b0726": 0.857,
            "b1241": 0.863,
            "b1276": 0.863,
            "b2935": 0.863,
            "b4025": 0.863,
        },
    }
    solution = double_gene_deletion(model, gene_list1=genes, processes=3)
    solution_one_process = double_gene_deletion(model, gene_list1=genes, processes=1)
    print(solution)
    for rxn_a, sub in iteritems(growth_dict):
        for rxn_b, growth in iteritems(sub):
            sol = solution.knockout[{rxn_a, rxn_b}]
            sol_one = solution_one_process.knockout[{rxn_a, rxn_b}]
            assert np.isclose(sol.growth, growth, atol=1e-3)
            assert np.isclose(sol_one.growth, growth, atol=1e-3)


# Double reaction deletion
def test_double_reaction_deletion_benchmark(large_model, benchmark):
    """Benchmark double reaction deletion."""
    reactions = large_model.reactions[1::100]
    benchmark(double_reaction_deletion, large_model, reaction_list1=reactions)


def test_double_reaction_deletion(model):
    """Test double reaction deletion."""
    reactions = ["FBA", "ATPS4r", "ENO", "FRUpts2"]
    growth_dict = {
        "FBA": {"ATPS4r": 0.135, "ENO": float("nan"), "FRUpts2": 0.704},
        "ATPS4r": {"ENO": float("nan"), "FRUpts2": 0.374},
        "ENO": {"FRUpts2": 0.0},
    }

    solution = double_reaction_deletion(model, reaction_list1=reactions, processes=3)
    solution_one_process = double_reaction_deletion(
        model, reaction_list1=reactions, processes=1
    )
    for (rxn_a, sub) in iteritems(growth_dict):
        for rxn_b, growth in iteritems(sub):
            sol = solution.knockout[{rxn_a, rxn_b}]
            sol_one = solution_one_process.knockout[{rxn_a, rxn_b}]
            if math.isnan(growth):
                assert math.isnan(sol.growth)
                assert math.isnan(sol_one.growth)
            else:
                assert np.isclose(sol.growth, growth, atol=1e-3)
                assert np.isclose(sol_one.growth, growth, atol=1e-3)
