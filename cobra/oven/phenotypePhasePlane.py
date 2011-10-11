from numpy import linspace, zeros, array, meshgrid, abs, empty, arange, int32
# not all plotting libraries may be installed
try:
    from cobra.external.ppmap import ppmap
except ImportError:
    ppmap = None
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
        self.segments = zeros(self.growth_rates.shape, dtype=int32)
        self.phases = []

    def plot(self):
        if mlab is not None:
            self.plot_mayavi()
        elif pyplot is not None:
            self.plot_matplotlib()
        else:
            raise (ImportError, "No suitable 3D plotting package found")

    def plot_matplotlib(self):
        """Use matplotlib to plot a phenotype phase plane in 3D.
        The resulting figure may be slow to interact with in real time,
        but will be easy to save as a vector figure.
        returns: maptlotlib 3d subplot object"""
        if pyplot is None:
            raise (ImportError, "Error importing matplotlib 3D plotting")
        # pick colors
        color_list = ('r', 'g', 'b', 'c', 'y', 'm', 'k')
        colors = empty(self.growth_rates.shape, dtype=str)
        n_segments = self.segments.max()
        if n_segments > 0 and n_segments <= len(color_list):
            for i in range(n_segments):
                colors[data.segments == (i + 1)] = color_list[i]
        else:
            colors[:, :] = 'b'
        figure = pyplot.figure()
        xgrid, ygrid = meshgrid(self.reaction1_fluxes, self.reaction2_fluxes)
        axes = figure.add_subplot(111, projection="3d")
        axes.plot_wireframe(xgrid, ygrid, self.growth_rates, color="black")
        axes.plot_surface(xgrid, ygrid, self.growth_rates, rstride=1,
            cstride=1, facecolors=colors, linewidth=0, antialiased=False)
        axes.set_xlabel(self.reaction1_name)
        axes.set_ylabel(self.reaction2_name)
        axes.set_zlabel("Growth rates")
        pyplot.show()
        return axes

    def plot_mayavi(self):
        """Use matplotlib to plot a phenotype phase plane in 3D.
        The resulting figure will be quick to interact with in real time,
        but might be difficult to save as a vector figure.
        returns: mlab figure object"""
        if mlab is None:
            raise (ImportError, "Error importing mayavi 3D plotting")
        figure = mlab.figure(bgcolor=(1, 1, 1))
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
            scalars=self.shadow_prices1 + self.shadow_prices2, resolution=1,
            representation="surface", opacity=0.75, figure=figure)
        # draw axes
        mlab.outline(extent=(0, max, 0, max, 0, max))
        mlab.axes(opacity=0, ranges=[0, xmax, 0, ymax, 0, zmax])
        mlab.xlabel(self.reaction1_name)
        mlab.ylabel(self.reaction2_name)
        mlab.zlabel("Growth rates")
        return figure

    def segment(self, threshold=0.005):
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
            i = self.segments.argmin() / self.segments.shape[0]
            j = self.segments.argmin() % self.segments.shape[0]
            # update the segment id for any point with a similar shadow price
            # to the current point
            self.segments[
                (abs(self.shadow_prices1 - self.shadow_prices1[i, j])
                    < threshold) *
                (abs(self.shadow_prices2 - self.shadow_prices2[i, j])
                    < threshold)] += segment_id
            # add the current point as one of the phases
            self.phases.append((
                (self.reaction1_fluxes[i], self.reaction2_fluxes[j]),
                self.shadow_prices1[i, j], self.shadow_prices2[i, j]))


def _calculate_subset(arguments):
    """Calculate a subset of the phenotype phase plane data.
    Store each result tuple as:
    (i, j, growth_rate, shadow_price1, shadow_price2)"""
    # unpack all arguments
    for key in arguments.iterkeys():
        exec("%s = arguments['%s']" % (key, key))
    results = []
    hot_start = None
    reaction1 = model.reactions[index1]
    reaction2 = model.reactions[index2]
    for a, flux1 in enumerate(reaction1_fluxes):
        i = i_list[a]
        reaction1.lower_bound = -1 * flux1 - tolerance
        reaction1.upper_bound = -1 * flux1 + tolerance
        for b, flux2 in enumerate(reaction2_fluxes):
            j = j_list[b]
            reaction2.lower_bound = -1 * flux2 - tolerance
            reaction2.upper_bound = -1 * flux2 + tolerance
            hot_start = model.optimize(the_problem=hot_start)
            if model.solution.status == "optimal":
                results.append((i, j, model.solution.f,
                    model.solution.y[index1],
                    model.solution.y[index2]))
            else:
                results.append((i, j, 0, 0, 0))
    return results


def calculate_phenotype_phase_plane(model,
        reaction1_name, reaction2_name,
        reaction1_range_max=20, reaction2_range_max=20,
        reaction1_npoints=50, reaction2_npoints=50,
        n_processes=1, tolerance=1e-6):
    """calculates the growth rates while varying the uptake rates for two
    reactions.

    returns: an object containing the growth rates for the uptake rates.
    To plot the result, call the plot function of the returned object.

    Example:
    data = calculate_phenotype_phase_plane(my_model, "EX_foo", "EX_bar")
    data.plot()
    """
    data = phenotypePhasePlaneData(
            str(reaction1_name), str(reaction2_name),
            reaction1_range_max, reaction2_range_max,
            reaction1_npoints, reaction2_npoints)
    # find the objects for the reactions and metabolites
    index1 = model.reactions.index(data.reaction1_name)
    index2 = model.reactions.index(data.reaction2_name)
    if n_processes > reaction1_npoints:
        n_processes = reaction1_npoints
    range_add = reaction1_npoints / n_processes
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
        arguments_list.append({"model": model.copy(),
            "index1": index1, "index2": index2,
            "reaction1_fluxes": r1_range,
            "reaction2_fluxes": data.reaction2_fluxes.copy(),
            "i_list": i_list, "j_list": j.copy(),
            "tolerance": tolerance})
    if n_processes > 1:
        results = list(ppmap(n_processes, _calculate_subset, arguments_list))
    else:
        results = [_calculate_subset(arguments_list[0])]
    for result_list in results:
        for result in result_list:
            i = result[0]
            j = result[1]
            data.growth_rates[i, j] = result[2]
            data.shadow_prices1[i, j] = result[3]
            data.shadow_prices2[i, j] = result[4]
    return data


if __name__ == "__main__":
    from time import time
    from os.path import join, split

    import cobra
    from cobra.io.sbml import create_cobra_model_from_sbml_file
    from cobra.test import data_directory

    n1 = 20
    n2 = 20
    the_file = join(data_directory, "salmonella.xml")
    model = create_cobra_model_from_sbml_file(the_file)
    print "sbml imported"
    start_time = time()
    data = calculate_phenotype_phase_plane(model, 'EX_glc__D_e',
        'EX_o2_e', reaction1_npoints=n1, reaction2_npoints=n2, n_processes=1)
    print "took %.2f seconds with 1 process" % (time() - start_time)
    start_time = time()
    data = calculate_phenotype_phase_plane(model, 'EX_glc__D_e',
        'EX_o2_e', reaction1_npoints=n1, reaction2_npoints=n2, n_processes=4)
    print "took %.2f seconds with 4 processes" % (time() - start_time)
    data.plot()
