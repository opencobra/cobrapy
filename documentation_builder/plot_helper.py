from matplotlib.pyplot import figure, xlim, ylim, gca, arrow, text, scatter
from mpl_toolkits.axes_grid.axislines import SubplotZero
from numpy import linspace, arange, sqrt, pi, sin, cos, sign
from IPython.display import set_matplotlib_formats

set_matplotlib_formats('png', 'pdf')


# axis style
def make_plot_ax():
    fig = figure(figsize=(6, 5))
    ax = SubplotZero(fig, 111)
    fig.add_subplot(ax)
    for direction in ["xzero", "yzero"]:
        ax.axis[direction].set_axisline_style("-|>")
        ax.axis[direction].set_visible(True)
    for direction in ["left", "right", "bottom", "top"]:
        ax.axis[direction].set_visible(False)
    xlim(-0.1, 2.1)
    ylim(xlim())
    ticks = [0.5 * i for i in range(1, 5)]
    labels = [str(i) if i == int(i) else "" for i in ticks]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.axis["yzero"].set_axis_direction("left")
    return ax


def plot_qp1():
    ax = make_plot_ax()
    ax.plot((0, 2), (2, 0), 'b')
    ax.plot([1], [1], 'bo')

    # circular grid
    for r in sqrt(2.) + 0.125 * arange(-11, 6):
        t = linspace(0., pi/2., 100)
        ax.plot(r * cos(t), r * sin(t), '-.', color="gray")


def plot_qp2():
    ax = make_plot_ax()
    ax.plot((0, 2), (2, 0), 'b')
    ax.plot([0.5], [1.5], 'bo')

    yrange = linspace(1, 2, 11)
    for r in (yrange ** 2 / 2. - yrange):
        t = linspace(-sqrt(2 * r + 1) + 0.000001,
                     sqrt(2 * r + 1) - 0.000001, 1000)
        ax.plot(abs(t), 1 + sqrt(2 * r + 1 - t ** 2) * sign(t), '-.',
                color="gray")


def plot_loop():
    figure(figsize=(10.5, 4.5), frameon=False)
    gca().axis("off")
    xlim(0.5, 3.5)
    ylim(0.7, 2.2)
    arrow_params = {"head_length": 0.08, "head_width": 0.1, "ec": "k",
                    "fc": "k"}
    text_params = {"fontsize": 25, "horizontalalignment": "center",
                   "verticalalignment": "center"}
    arrow(0.5, 1, 0.85, 0, **arrow_params)  # EX_A
    arrow(1.5, 1, 0.425, 0.736, **arrow_params)  # v1
    arrow(2.04, 1.82, 0.42, -0.72, **arrow_params)  # v2
    arrow(2.4, 1, -0.75, 0, **arrow_params)  # v3
    arrow(2.6, 1, 0.75, 0, **arrow_params)
    # reaction labels
    text(0.9, 1.15, "EX_A", **text_params)
    text(1.6, 1.5, r"v$_1$", **text_params)
    text(2.4, 1.5, r"v$_2$", **text_params)
    text(2, 0.85, r"v$_3$", **text_params)
    text(2.9, 1.15, "DM_C", **text_params)
    # metabolite labels
    scatter(1.5, 1, s=250, color='#c994c7')
    text(1.5, 0.9, "A", **text_params)
    scatter(2, 1.84, s=250, color='#c994c7')
    text(2, 1.95, "B", **text_params)
    scatter(2.5, 1, s=250, color='#c994c7')
    text(2.5, 0.9, "C", **text_params)
