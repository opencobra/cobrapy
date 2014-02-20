import cobra.io
import cobra.test
from optGpSampler.sampler import Sampler
core_model = cobra.io.load_matlab_model(cobra.test.ecoli_core_mat)
core_model.to_array_based_model()
sampler = Sampler(core_model.S.tolil(), core_model.lower_bounds, core_model.upper_bounds)
sampler.setSolver("cplex")
sampler.nSamples = 10000
# As of 2/19/14 there is a bug where optGpSampler leaves out the last reaction.
# This next line will fix it.
sampler.setReactionIdsToKeep("all")
samples = sampler.sample()
