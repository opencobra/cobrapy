# -*- coding: utf-8 -*-

import locale

import cobra
import cobra.test


if __name__ == "__main__":
    locale.setlocale(locale.LC_ALL, "")
    model = cobra.test.create_test_model("ecoli")
    solution = cobra.flux_analysis.double_reaction_deletion(model,
            reaction_list1=model.reactions[:200], num_cpu=4)
