# This example demonstrates reading in a COBRA SBML model xml file
# performing a simple optimization, and then saving the model as an SBML xml file.
#
# The model is from the Salmonella enterica Typhimurium LT2 publication in
# 2011 in BMC Sys Bio 5:8
#
# The simulated growth rate should be ~0.478.
#
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file
from cobra.test import salmonella_sbml
sbml_out_file = 'salmonella.out.xml'
sbml_level = 2
sbml_version = 1  # Writing version 4 is not completely supported.
#Read in the sbml file.
cobra_model = create_cobra_model_from_sbml_file(salmonella_sbml)
# Run the optimization for the objective reaction and medium composition
# set in the file.
cobra_model.optimize(solver='glpk')  # 'gurobi' or 'cplex' are also supported
print('\nSimulated growth rate is %1.3f' % cobra_model.solution.f)

#Save the model to an SBML file
print('\nConverting Model to SBML and saving as ' + sbml_out_file)
write_cobra_model_to_sbml_file(cobra_model, sbml_out_file, sbml_level,
                               sbml_version)
