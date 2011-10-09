from numpy import linspace, zeros, array, meshgrid, abs
# not all plotting libraries may be installed
try:
    from mayavi import mlab
except ImportError:
    mlab = None
try:
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import axes3d
except ImportError:
    pyplot = None
    axes3d = None


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
            reaction1_npoints)
        self.growth_rates = zeros((reaction1_npoints, reaction2_npoints))
        self.shadow_prices1 = zeros((reaction1_npoints, reaction2_npoints))
        self.shadow_prices2 = zeros((reaction1_npoints, reaction2_npoints))

    def plot(self):
        self.plot_matplotlib()
        # self.plot_mayavi()

    def plot_matplotlib(self):
        if pyplot is None:
            raise (ImportError, "Error importing matplotlib 3D plotting")
        figure = pyplot.figure()
        xgrid, ygrid = meshgrid(self.reaction1_fluxes, self.reaction2_fluxes)
        axes = figure.add_subplot(111, projection="3d")
        axes.plot_wireframe(xgrid, ygrid, self.growth_rates, color="black")
        axes.plot_surface(xgrid, ygrid, self.growth_rates, alpha=0.75)
        axes.set_xlabel(self.reaction1_name)
        axes.set_ylabel(self.reaction2_name)
        axes.set_zlabel("Growth rates")
        pyplot.show()

    def plot_mayavi(self):
        if mlab is None:
            raise (ImportError, "Error importing mayavi 3D plotting")
        mlab.figure(bgcolor=(1, 1, 1))
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
            representation="wireframe", color=(0, 0, 0))
        mlab.surf(xgrid * xscale, ygrid * yscale, self.growth_rates * zscale,
            representation="surface", opacity=0.75)
        # draw axes
        mlab.outline(extent=(0, max, 0, max, 0, max))
        mlab.axes(opacity=0, ranges=[0, xmax, 0, ymax, 0, zmax])
        mlab.xlabel(self.reaction1_name)
        mlab.ylabel(self.reaction2_name)
        mlab.zlabel("Growth rates")


# TODO - parallelize
def calculate_phenotype_phase_plane(model,
        reaction1_name, reaction2_name,
        reaction1_range_max=20, reaction2_range_max=20,
        reaction1_npoints=50, reaction2_npoints=50,
        tolerance=1e-6):
    """calculates the growth rates while varying the uptake rates for two
    reactions.

    returns: an object containing the growth rates for the uptake rates.
    To plot the result, call the plot function of the returned object.

    Example:
    data = calculate_phenotype_phase_plane(my_model, "EX_foo", "EX_bar")
    data.plot()
    """
    data = phenotypePhasePlaneData(str(reaction1_name), str(reaction2_name),
            reaction1_range_max, reaction2_range_max,
            reaction1_npoints, reaction2_npoints)
    # find the objects for the reactions and metabolites
    reaction1 = model.reactions[model.reactions.index(data.reaction1_name)]
    reaction2 = model.reactions[model.reactions.index(data.reaction2_name)]
    metabolite1 = model.metabolites[model.metabolites.index(
        data.reaction1_name[3:])]
    metabolite2 = model.metabolites[model.metabolites.index(
        data.reaction2_name[3:])]
    # create a hot_start to help make solving faster
    hot_start = None
    for i in range(data.reaction1_npoints):
        reaction1.lower_bound = -1 * data.reaction1_fluxes[i] - tolerance
        reaction1.upper_bound = -1 * data.reaction1_fluxes[i] + tolerance
        for j in range(data.reaction2_npoints):
            reaction2.lower_bound = -1 * data.reaction1_fluxes[j] - tolerance
            reaction2.upper_bound = -1 * data.reaction1_fluxes[j] + tolerance
            hot_start = model.optimize(the_problem=hot_start)
            if model.solution.status == "optimal":
                data.growth_rates[i, j] = model.solution.f
                data.shadow_prices1[i, j] =
                    abs(model.solution.y_dict[data.reaction1_name])
                data.shadow_prices2[i, j] =
                    abs(model.solution.y_dict[data.reaction2_name])
            else:
                data.growth_rates[i, j] = 0
                data.shadow_prices1[i, j] = 0
                data.shadow_prices2[i, j] = 0
    return data
if __name__ == "__main__":
    import cobra
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    model = create_cobra_model_from_sbml_file(
        r"C:\Python27\Lib\site-packages\cobra\examples\files\salmonella.xml")
    data = calculate_phenotype_phase_plane(model, 'EX_glc__D_e', 'EX_o2_e',
        reaction1_npoints=10, reaction2_npoints=10)
    data.plot()
