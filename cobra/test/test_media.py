# -*- coding: utf-8 -*-

from __future__ import absolute_import
import pandas as pd
import cobra.media.media as media


class TestExchangeDetection:

    def test_external_compartment(self, model):
        assert media.external_compartment(model) == "e"

    def test_exchanges(self, model):
        ex = media.exchanges(model)
        assert all(r.id.startswith("EX_") for r in ex)


class TestMinimalMedia:

    def test_medium_linear(self, model):
        medium = media.minimal_medium(model, 0.8)
        assert len(medium) <= 4
        assert all(medium > 1e-6)

    def test_medium_mip(self, model):
        medium = media.minimal_medium(model, 0.8, minimize_components=True)
        assert len(medium) <= 4
        assert all(medium > 1e-6)

        # Anaerobic growth
        medium = media.minimal_medium(model, 0.1, minimize_components=True)
        assert len(medium) <= 3
        assert all(medium > 1e-6)

    def test_medium_alternative_mip(self, model):
        medium = media.minimal_medium(model, 0.8, minimize_components=5,
                                      open_exchanges=True)
        assert isinstance(medium, pd.DataFrame)
        assert medium.shape[1] == 5
        assert all((medium > 0).sum() <= 4)
        assert all(medium.sum(axis=1) > 1e-6)

    def test_benchmark_medium_linear(self, model, benchmark):
        benchmark(media.minimal_medium, model, 0.8)

    def test_benchmark_medium_mip(self, model, benchmark):
        benchmark(media.minimal_medium, model, 0.8, True)

    def test_medium_exports(self, model):
        medium = media.minimal_medium(model, 0.8, exports=True,
                                      minimize_components=True)
        assert len(medium) > 4
        assert any(medium < -1e-6)

    def test_open_exchanges(self, model):
        model.reactions.EX_glc__D_e.bounds = 0, 0
        medium = media.minimal_medium(model, 0.8)
        assert medium is None
        medium = media.minimal_medium(model, 0.8, minimize_components=True)
        assert medium is None

        medium = media.minimal_medium(model, 0.8, open_exchanges=True)
        assert len(medium) >= 3
        medium = media.minimal_medium(model, 0.8, open_exchanges=100)
        assert len(medium) >= 3
