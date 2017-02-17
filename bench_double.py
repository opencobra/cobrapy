# -*- coding: utf-8 -*-

import cobra
import cobra.test


if __name__ == "__main__":
    model = cobra.test.create_test_model("textbook")
    solution = cobra.flux_analysis.double_reaction_deletion(model, num_cpu=2)
