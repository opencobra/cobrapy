"""Test functionalities of gap filling."""

from cobra.core import Metabolite, Model, Reaction
from cobra.flux_analysis.gapfilling import GapFiller, gapfill


def test_gapfilling(salmonella: Model) -> None:
    """Test Gapfilling."""
    m = Model()
    m.add_metabolites([Metabolite(m_id) for m_id in ["a", "b", "c"]])
    exa = Reaction("EX_a")
    exa.add_metabolites({m.metabolites.a: 1})
    b2c = Reaction("b2c")
    b2c.add_metabolites({m.metabolites.b: -1, m.metabolites.c: 1})
    dmc = Reaction("DM_c")
    dmc.add_metabolites({m.metabolites.c: -1})
    m.add_reactions([exa, b2c, dmc])
    m.objective = "DM_c"

    universal = Model()
    a2b = Reaction("a2b")
    a2d = Reaction("a2d")
    universal.add_reactions([a2b, a2d])
    a2b.build_reaction_from_string("a --> b", verbose=False)
    a2d.build_reaction_from_string("a --> d", verbose=False)

    # # GrowMatch
    # result = gapfilling.growMatch(m, universal)[0]
    result = gapfill(m, universal)[0]
    assert len(result) == 1
    assert result[0].id == "a2b"

    # # SMILEY
    # result = gapfilling.SMILEY(m, "b", universal)[0]
    with m:
        m.objective = m.add_boundary(m.metabolites.b, type="demand")
        result = gapfill(m, universal)[0]
        assert len(result) == 1
        assert result[0].id == "a2b"

    # # 2 rounds of GrowMatch with exchange reactions
    # result = gapfilling.growMatch(m, None, ex_rxns=True, iterations=2)
    result = gapfill(m, None, exchange_reactions=True, iterations=2)
    assert len(result) == 2
    assert len(result[0]) == 1
    assert len(result[1]) == 1
    assert {i[0].id for i in result} == {"EX_b", "EX_c"}

    # # Gapfilling solution adds metabolites not present in original model
    # test for when demand = T
    # a demand reaction must be added to clear new metabolite
    universal_noDM = Model()
    a2b = Reaction("a2b")
    universal_noDM.add_reactions([a2b])
    a2b.build_reaction_from_string("a --> b + d", verbose=False)
    result = gapfill(
        m, universal_noDM, exchange_reactions=False, demand_reactions=True
    )[0]
    # add reaction a2b and demand reaction to clear met d
    assert len(result) == 2
    assert "a2b" in [x.id for x in result]

    # test for when demand = False
    # test for when metabolites are added to the model and
    # must be cleared by other reactions in universal model
    # (i.e. not necessarily a demand reaction)
    universal_withDM = universal_noDM.copy()
    d_dm = Reaction("d_dm")
    universal_withDM.add_reactions([d_dm])
    d_dm.build_reaction_from_string("d -->", verbose=False)
    result = gapfill(
        m, universal_withDM, exchange_reactions=False, demand_reactions=False
    )[0]
    assert len(result) == 2
    assert "a2b" in [x.id for x in result]

    # somewhat bigger model
    universal = Model("universal_reactions")
    with salmonella as model:
        for i in [i.id for i in model.metabolites.f6p_c.reactions]:
            reaction = model.reactions.get_by_id(i)
            universal.add_reactions([reaction.copy()])
            model.remove_reactions([reaction])
        gf = GapFiller(
            model, universal, penalties={"TKT2": 1e3}, demand_reactions=False
        )
        solution = gf.fill()
        assert "TKT2" not in {r.id for r in solution[0]}
        assert gf.validate(solution[0])
