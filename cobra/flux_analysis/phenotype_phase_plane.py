# -*- coding: utf-8 -*-

from __future__ import absolute_import, division

from warnings import warn
from collections import defaultdict
from operator import itemgetter
from itertools import product
from six import iteritems
from multiprocessing import Pool
import pandas as pd
from optlang.interface import OPTIMAL

from numpy import (
    nan, abs, arange, dtype, empty, int32, linspace, meshgrid, unravel_index,
    zeros, array)

from cobra.solvers import get_solver_name, solver_dict
import cobra.util.solver as sutil
from cobra.flux_analysis import flux_variability_analysis as fva
from cobra.exceptions import OptimizationError

# attempt to import plotting libraries
try:
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import axes3d
except (ImportError, RuntimeError):
    pyplot = None
    axes3d = None
mlab = None  # mayavi may crash python
try:  # for prettier colors
    from palettable.colorbrewer import get_map
except (ImportError, RuntimeError):
    try:
        from brewer2mpl import get_map
    except ImportError:
        get_map = None


class phenotypePhasePlaneData:
    """class to hold results of a phenotype phase plane analysis"""

    def __init__(self,
                 reaction1_name, reaction2_name,
                 reaction1_range_max, reaction2_range_max,
                 reaction1_npoints, reaction2_npoints):
        self.reaction1_name = reaction1_name
        self.reaction2_name = reaction2_name
        self.reaction1_range_max = reaction1_range_max
        self.reaction2_range_max = reaction2_range_max
        self.reaction1_npoints = reaction1_npoints
        self.reaction2_npoints = reaction2_npoints
        self.reaction1_fluxes = linspace(0, reaction1_range_max,
                                         reaction1_npoints)
        self.reaction2_fluxes = linspace(0, reaction2_range_max,
                                         reaction2_npoints)
        self.growth_rates = zeros((reaction1_npoints, reaction2_npoints))
        self.shadow_prices1 = zeros((reaction1_npoints, reaction2_npoints))
        self.shadow_prices2 = zeros((reaction1_npoints, reaction2_npoints))
        self.segments = zeros(self.growth_rates.shape, dtype=int32)
        self.phases = []

    def plot(self):
        """plot the phenotype phase plane in 3D using any available backend"""
        if pyplot is not None:
            self.plot_matplotlib()
        elif mlab is not None:
            self.plot_mayavi()
        else:
            raise ImportError("No suitable 3D plotting package found")

    def plot_matplotlib(self, theme="Paired", scale_grid=False):
        """Use matplotlib to plot a phenotype phase plane in 3D.

        theme: color theme to use (requires palettable)

        returns: maptlotlib 3d subplot object"""
        if pyplot is None:
            raise ImportError("Error importing matplotlib 3D plotting")
        colors = empty(self.growth_rates.shape, dtype=dtype((str, 7)))
        n_segments = self.segments.max()
        # pick colors
        color_list = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
                      '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00',
                      '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928']
        if get_map is not None:
            try:
                color_list = get_map(theme, 'Qualitative',
                                     n_segments).hex_colors
            except ValueError:
                from warnings import warn
                warn('palettable could not be used for this number of phases')
        if n_segments > len(color_list):
            from warnings import warn
            warn("not enough colors to color all detected phases")
        if n_segments > 0 and n_segments <= len(color_list):
            for i in range(n_segments):
                colors[self.segments == (i + 1)] = color_list[i]
        else:
            colors[:, :] = 'b'
        if scale_grid:
            # grid wires should not have more than ~20 points
            xgrid_scale = int(self.reaction1_npoints / 20)
            ygrid_scale = int(self.reaction2_npoints / 20)
        else:
            xgrid_scale, ygrid_scale = (1, 1)
        figure = pyplot.figure()
        xgrid, ygrid = meshgrid(self.reaction1_fluxes, self.reaction2_fluxes)
        axes = figure.add_subplot(111, projection="3d")
        xgrid = xgrid.transpose()
        ygrid = ygrid.transpose()
        axes.plot_surface(xgrid, ygrid, self.growth_rates, rstride=1,
                          cstride=1, facecolors=colors, linewidth=0,
                          antialiased=False)
        axes.plot_wireframe(xgrid, ygrid, self.growth_rates, color="black",
                            rstride=xgrid_scale, cstride=ygrid_scale)
        axes.set_xlabel(self.reaction1_name, size="x-large")
        axes.set_ylabel(self.reaction2_name, size="x-large")
        axes.set_zlabel("Growth rate", size="x-large")
        axes.view_init(elev=30, azim=-135)
        figure.set_tight_layout(True)
        return axes

    def plot_mayavi(self):
        """Use mayavi to plot a phenotype phase plane in 3D.
        The resulting figure will be quick to interact with in real time,
        but might be difficult to save as a vector figure.
        returns: mlab figure object"""
        from mayavi import mlab
        figure = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
        figure.name = "Phenotype Phase Plane"
        max = 10.0
        xmax = self.reaction1_fluxes.max()
        ymax = self.reaction2_fluxes.max()
        zmax = self.growth_rates.max()
        xgrid, ygrid = meshgrid(self.reaction1_fluxes, self.reaction2_fluxes)
        xgrid = xgrid.transpose()
        ygrid = ygrid.transpose()
        xscale = max / xmax
        yscale = max / ymax
        zscale = max / zmax
        mlab.surf(xgrid * xscale, ygrid * yscale, self.growth_rates * zscale,
                  representation="wireframe", color=(0, 0, 0), figure=figure)
        mlab.mesh(xgrid * xscale, ygrid * yscale, self.growth_rates * zscale,
                  scalars=self.shadow_prices1 + self.shadow_prices2,
                  resolution=1, representation="surface", opacity=0.75,
                  figure=figure)
        # draw axes
        mlab.outline(extent=(0, max, 0, max, 0, max))
        mlab.axes(opacity=0, ranges=[0, xmax, 0, ymax, 0, zmax])
        mlab.xlabel(self.reaction1_name)
        mlab.ylabel(self.reaction2_name)
        mlab.zlabel("Growth rates")
        return figure

    def segment(self, threshold=0.01):
        """attempt to segment the data and identify the various phases"""
        self.segments *= 0
        # each entry in phases will consist of the following tuple
        # ((x, y), shadow_price1, shadow_price2)
        self.phases = []
        # initialize the area to be all False
        covered_area = (self.growth_rates * 0 == 1)
        # as long as part of the area has not been covered
        segment_id = 0
        while self.segments.min() == 0:
            segment_id += 1
            # i and j are indices for a current point which has not been
            # assigned a segment yet

            i, j = unravel_index(self.segments.argmin(), self.segments.shape)
            # update the segment id for any point with a similar shadow price
            # to the current point
            d1 = abs(self.shadow_prices1 - self.shadow_prices1[i, j])
            d2 = abs(self.shadow_prices2 - self.shadow_prices2[i, j])
            self.segments[(d1 < threshold) * (d2 < threshold)] += segment_id
            # add the current point as one of the phases
            self.phases.append((
                (self.reaction1_fluxes[i], self.reaction2_fluxes[j]),
                self.shadow_prices1[i, j], self.shadow_prices2[i, j]))


def _calculate_subset(arguments):
    """Calculate a subset of the phenotype phase plane data.
    Store each result tuple as:
    (i, j, growth_rate, shadow_price1, shadow_price2)"""

    model = arguments["model"]
    reaction1_fluxes = arguments["reaction1_fluxes"]
    reaction2_fluxes = arguments["reaction2_fluxes"]
    metabolite1_name = arguments["metabolite1_name"]
    metabolite2_name = arguments["metabolite2_name"]
    index1 = arguments["index1"]
    index2 = arguments["index2"]
    i_list = arguments["i_list"]
    j_list = arguments["j_list"]
    tolerance = arguments["tolerance"]
    solver = solver_dict[arguments["solver"]]

    results = []
    reaction1 = model.reactions[index1]
    reaction2 = model.reactions[index2]
    problem = solver.create_problem(model)
    solver.solve_problem(problem)
    for a, flux1 in enumerate(reaction1_fluxes):
        i = i_list[a]
        # flux is actually negative for uptake. Also some solvers require
        # float instead of numpy.float64
        flux1 = float(-1 * flux1)
        # change bounds on reaction 1
        solver.change_variable_bounds(problem, index1, flux1 - tolerance,
                                      flux1 + tolerance)
        for b, flux2 in enumerate(reaction2_fluxes):
            j = j_list[b]
            flux2 = float(-1 * flux2)  # same story as flux1
            # change bounds on reaction 2
            solver.change_variable_bounds(problem, index2, flux2 - tolerance,
                                          flux2 + tolerance)
            # solve the problem and save results
            solver.solve_problem(problem)
            solution = solver.format_solution(problem, model)
            if solution is not None and solution.status == "optimal":
                results.append((i, j, solution.f,
                                solution.y_dict[metabolite1_name],
                                solution.y_dict[metabolite2_name]))
            else:
                results.append((i, j, 0, 0, 0))
            # reset reaction 2 bounds
            solver.change_variable_bounds(problem, index2,
                                          float(reaction2.lower_bound),
                                          float(reaction2.upper_bound))
        # reset reaction 1 bounds
        solver.change_variable_bounds(problem, index1,
                                      float(reaction1.lower_bound),
                                      float(reaction1.upper_bound))
    return results


def calculate_phenotype_phase_plane(
        model, reaction1_name, reaction2_name,
        reaction1_range_max=20, reaction2_range_max=20,
        reaction1_npoints=50, reaction2_npoints=50,
        solver=None, n_processes=1, tolerance=1e-6):
    """calculates the growth rates while varying the uptake rates for two
    reactions.

    :returns: a `phenotypePhasePlaneData` object containing the growth rates
    for the uptake rates. To plot the
    result, call the plot function of the returned object.

    :Example:
    >>> import cobra.test
    >>> model = cobra.test.create_test_model("textbook")
    >>> ppp = calculate_phenotype_phase_plane(model, "EX_glc__D_e", "EX_o2_e")
    >>> ppp.plot()
    """
    warn('calculate_phenotype_phase_plane is deprecated, consider using '
         'production_envelope instead', DeprecationWarning)
    if solver is None:
        solver = get_solver_name()
    data = phenotypePhasePlaneData(
        str(reaction1_name), str(reaction2_name),
        reaction1_range_max, reaction2_range_max,
        reaction1_npoints, reaction2_npoints)
    # find the objects for the reactions and metabolites
    index1 = model.reactions.index(data.reaction1_name)
    index2 = model.reactions.index(data.reaction2_name)
    metabolite1_name = list(model.reactions[index1]._metabolites)[0].id
    metabolite2_name = list(model.reactions[index2]._metabolites)[0].id
    if n_processes > reaction1_npoints:  # limit the number of processes
        n_processes = reaction1_npoints
    range_add = reaction1_npoints // n_processes
    # prepare the list of arguments for each _calculate_subset call
    arguments_list = []
    i = arange(reaction1_npoints)
    j = arange(reaction2_npoints)
    for n in range(n_processes):
        start = n * range_add
        if n != n_processes - 1:
            r1_range = data.reaction1_fluxes[start:start + range_add]
            i_list = i[start:start + range_add]
        else:
            r1_range = data.reaction1_fluxes[start:]
            i_list = i[start:]
        arguments_list.append({
            "model": model,
            "index1": index1, "index2": index2,
            "metabolite1_name": metabolite1_name,
            "metabolite2_name": metabolite2_name,
            "reaction1_fluxes": r1_range,
            "reaction2_fluxes": data.reaction2_fluxes.copy(),
            "i_list": i_list, "j_list": j.copy(),
            "tolerance": tolerance, "solver": solver})
    if n_processes > 1:
        p = Pool(n_processes)
        results = list(p.map(_calculate_subset, arguments_list))
    else:
        results = [_calculate_subset(arguments_list[0])]
    for result_list in results:
        for result in result_list:
            i = result[0]
            j = result[1]
            data.growth_rates[i, j] = result[2]
            data.shadow_prices1[i, j] = result[3]
            data.shadow_prices2[i, j] = result[4]
    data.segment()
    return data


def production_envelope(model, reactions, objective=None, c_source=None,
                        points=20, solver=None):
    """Calculate the objective value conditioned on all combinations of
    fluxes for a set of chosen reactions

    The production envelope can be used to analyze a models ability to
    produce a given compound conditional on the fluxes for another set of
    reaction, such as the uptake rates. The model is alternately optimize
    with respect to minimizing and maximizing the objective and record the
    obtained fluxes. Ranges to compute production is set to the effective
    bounds, i.e. the minimum / maximum fluxes that can be obtained given
    current reaction bounds.

    Parameters
    ----------
    model : cobra.Model
        The model to compute the production envelope for.
    reactions : list or string
        A list of reactions, reaction identifiers or single reaction
    objective : string, dict, model.solver.interface.Objective
        The objective (reaction) to use for the production envelope. Use the
        model's current objective is left missing.
    c_source : cobra.Reaction or string
       A reaction or reaction identifier that is the source of carbon for
       computing carbon (mol carbon in output over mol carbon in input) and
       mass yield (gram product over gram output). Only objectives with a
       carbon containing input and output metabolite is supported.
    points : int
       The number of points to calculate production for.
    solver : string
       The solver to use - only here for consistency with older
       implementations (this argument will be removed in the future). The
       solver should be set using `model.solver` directly. Only optlang
       based solvers are supported.

    Returns
    -------
    pandas.DataFrame
        A data frame with one row per evaluated point and

        - reaction id : one column per input reaction indicating the flux at
          each given point,
        - flux: the objective flux

        - carbon_yield: if carbon source is defined and the product is a
          single metabolite (mol carbon product per mol carbon feeding source)

        - mass_yield: if carbon source is defined and the product is a
          single metabolite (gram product per 1 g of feeding source)

        - direction: the direction of the optimization.

        Only points that give a valid solution are returned.

    Examples
    --------
    >>> import cobra.test
    >>> from cobra.flux_analysis import production_envelope
    >>> model = cobra.test.create_test_model("textbook")
    >>> production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
    """
    legacy, solver = sutil.choose_solver(model, solver)
    if legacy:
        raise ValueError('production_envelope is only implemented for optlang '
                         'based solver interfaces.')
    reactions = model.reactions.get_by_any(reactions)
    objective = model.solver.objective if objective is None else objective
    result = None
    with model:
        model.objective = objective
        if c_source is None:
            c_input = get_c_input(model)
        else:
            c_input = model.reactions.get_by_any(c_source)[0]
        objective_reactions = list(sutil.linear_reaction_coefficients(model))
        if len(objective_reactions) != 1:
            raise ValueError('cannot calculate yields for objectives with '
                             'multiple reactions')
        carbon_io = c_input, objective_reactions[0]
        min_max = fva(model, reactions, fraction_of_optimum=0)
        grid = [linspace(min_max.minimum[rxn.id], min_max.maximum[rxn.id],
                         points, endpoint=True) for rxn in reactions]
        grid_list = list(product(*grid))
        result = envelope_for_points(model, reactions, grid_list, carbon_io)

    return pd.DataFrame(result)


def envelope_for_points(model, reactions, grid, carbon_io):
    results = defaultdict(list)
    for direction in ('minimum', 'maximum'):
        sense = "min" if direction == "minimum" else "max"
        for point in grid:
            with model:
                model.solver.objective.direction = sense
                for reaction, coordinate in zip(reactions, point):
                    reaction.bounds = (coordinate, coordinate)
                model.solver.optimize()
                if model.solver.status == OPTIMAL:
                    for reaction, coordinate in zip(reactions, point):
                        results[reaction.id].append(coordinate)
                    results['direction'].append(direction)
                    results['flux'].append(model.solver.objective.value)
                    if carbon_io[0] is not None:
                        results['carbon_yield'].append(carbon_yield(carbon_io))
                        results['mass_yield'].append(mass_yield(carbon_io))
    for key, value in results.items():
        results[key] = array(value)
    if carbon_io[0] is not None:
        results['carbon_source'] = carbon_io[0].id
    return results


def carbon_yield(c_input_output):
    """ mol product per mol carbon input

    Returns
    -------
    float
        the mol carbon atoms in the product (as defined by the model
        objective) divided by the mol carbon in the input reactions (as
        defined by the model medium) or zero in case of division by zero
        arises
    """

    c_input, c_output = c_input_output
    if c_input is None:
        return nan
    carbon_input_flux = total_carbon_flux(c_input, consumption=True)
    carbon_output_flux = total_carbon_flux(c_output, consumption=False)
    try:
        return carbon_output_flux / carbon_input_flux
    except ZeroDivisionError:
        return nan


def mass_yield(c_input_output):
    """Gram product divided by gram of carbon input source

    Parameters
    ----------
    c_input_output : tuple
        Two reactions, the one that feeds carbon to the system and the one
        that produces carbon containing compound.

    Returns
    -------
    float
        gram product per 1 g of feeding source
    """
    c_input, c_output = c_input_output
    if input is None:
        return nan
    try:
        c_source, source_flux = single_flux(c_input, consumption=True)
        c_product, product_flux = single_flux(c_output, consumption=False)
    except ValueError:
        return nan
    mol_prod_mol_src = product_flux / source_flux
    x = mol_prod_mol_src * c_product.formula_weight
    return x / c_source.formula_weight


def total_carbon_flux(reaction, consumption=True):
    """summed product carbon flux for a reaction

    Parameters
    ----------
    reaction : Reaction
        the reaction to carbon return flux for
    consumption : bool
        flux for consumed metabolite, else produced

    Returns
    -------
    float
        reaction flux multiplied by number of carbon for the products of the
        reaction
    """
    direction = 1 if consumption else -1
    c_flux = [reaction.flux * coeff * met.elements.get('C', 0) * direction
              for met, coeff in reaction.metabolites.items()]
    return sum([flux for flux in c_flux if flux > 0])


def single_flux(reaction, consumption=True):
    """flux into single product for a reaction

    only defined for reactions with single products

    Parameters
    ----------
    reaction : Reaction
        the reaction to product flux for
    consumption : bool
        flux for consumed metabolite, else produced

    Returns
    -------
    tuple
        metabolite, flux for the metabolite
    """
    if len(list(reaction.metabolites)) != 1:
        raise ValueError('product flux only defined for single metabolite '
                         'reactions')
    met, coeff = next(iteritems(reaction.metabolites))
    direction = 1 if consumption else -1
    return met, reaction.flux * coeff * direction


def get_c_input(model):
    """ carbon source reactions

    Returns
    -------
    Reaction
       The medium reaction with highest input carbon flux
    """
    try:
        model.solver.optimize()
        sutil.assert_optimal(model)
    except OptimizationError:
        return None

    reactions = model.reactions.get_by_any(list(model.medium))
    reactions_fluxes = [(rxn, total_carbon_flux(rxn, consumption=True))
                        for rxn in reactions]
    source_reactions = [(rxn, c_flux) for rxn, c_flux
                        in reactions_fluxes if c_flux > 0]
    try:
        return max(source_reactions, key=itemgetter(1))[0]
    except ValueError:
        return None
