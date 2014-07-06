
## Escher and cobrapy

# cobrapy integrates well with the
# [escher](https://github.com/zakandrewking/visbio) package. If it has been
# installed, escher makes downloading and viewing metabolic maps from a
# [repository](https://github.com/zakandrewking/escher/tree/gh-pages/maps)
# extremely simple. For more information, view the escher documentation.
# 
# This example is also available as an IPython [notebook](http://nbviewer.ipyth
# on.org/github/opencobra/cobrapy/blob/master/documentation_builder/escher.ipyn
# b)

from escher import Builder
import cobra.test

model = cobra.test.create_test_model(cobra.test.ecoli_pickle)
cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
wt_solution = model.solution.x_dict

# mutant flux has PGI knocked out
model.reactions.PGI.lower_bound = 0
model.reactions.PGI.upper_bound = 0
cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)
mutant_solution = model.solution.x_dict

wt_map = Builder("iJO1366_central_metabolism", reaction_data=wt_solution)
wt_map.display_in_notebook()

wt_map = Builder("iJO1366_central_metabolism", reaction_data=mutant_solution)
wt_map.display_in_notebook()

# In a non-notebook environment, metabolic maps can still be viewed by running

                mutant_map.display_in_browser()
                