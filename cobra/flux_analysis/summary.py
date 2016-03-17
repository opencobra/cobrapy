from six.moves import zip_longest
from six import iterkeys, print_, text_type

import pandas as pd

from .variability import flux_variability_analysis


def format_long_string(string, max_length):
    if len(string) > max_length:
        string = string[:max_length - 3]
        string += '...'
    return string


def metabolite_summary(met, threshold=0.01, fva=False, **solver_args):
    """Print a summary of the reactions which produce and consume this
    metabolite

    threshold: float
    a value below which to ignore reaction fluxes

    fva: float (0->1), or None
    Whether or not to include flux variability analysis in the output.
    If given, fva should be a float between 0 and 1, representing the
    fraction of the optimum objective to be searched.

    """

    def rxn_summary(r):
        return {
            'id': r.id,
            'flux': r.x * r.metabolites[met],
            'reaction': r.reaction,
        }

    flux_summary = pd.DataFrame((rxn_summary(r) for r in met.reactions))
    assert flux_summary.flux.sum() < 1E-6, "Error in flux balance"
    producing = flux_summary[flux_summary.flux > 0].copy()
    consuming = flux_summary[flux_summary.flux < 0].copy()

    for df in [producing, consuming]:
        df['percent'] = df.flux / df.flux.sum()
        df.drop(df[df['percent'] < threshold].index, axis=0,
                inplace=True)

    producing.sort_values('percent', ascending=False, inplace=True)
    consuming.sort_values('percent', ascending=False, inplace=True)

    if not fva:

        producing.flux = producing.flux.apply(
            lambda x: '{:6.2g}'.format(x))
        consuming.flux = consuming.flux.apply(
            lambda x: '{:6.2g}'.format(x))

        flux_len = 6

    else:
        fva_results = pd.DataFrame(
            flux_variability_analysis(met.model, met.reactions,
                                      fraction_of_optimum=fva,
                                      **solver_args)).T
        half_span = (fva_results.maximum - fva_results.minimum) / 2
        median = fva_results.minimum + half_span

        producing.flux = producing.id.apply(
            lambda x, median=median, err=half_span:
            u'{0:0.2f} \u00B1 {1:0.2f}'.format(median[x], err[x]))
        consuming.flux = consuming.id.apply(
            lambda x, median=median, err=half_span:
            u'{0:0.2f} \u00B1 {1:0.2f}'.format(median[x], err[x]))

        flux_len = max(producing.flux.apply(len).max(),
                       consuming.flux.apply(len).max()) + 1

    for df in [producing, consuming]:

        df['reaction'] = df['reaction'].map(
            lambda x: format_long_string(x, 52))
        df['id'] = df['id'].map(
            lambda x: format_long_string(x, 8))

    head = "PRODUCING REACTIONS -- " + format_long_string(met.name, 55)
    print_(head)
    print_("-" * len(head))
    print_(("{0:^6} {1:>" + str(flux_len) + "} {2:>8} {3:^54}").format(
        '%', 'FLUX', 'RXN ID', 'REACTION'))

    for row in producing.iterrows():
        print_((u"{0.percent:6.1%} {0.flux:>" + str(flux_len) +
                "} {0.id:>8} {0.reaction:>54}").format(row[1]))

    print_()
    print_("CONSUMING REACTIONS -- " + format_long_string(met.name, 55))
    print_("-" * len(head))
    print_(("{0:^6} {1:>" + str(flux_len) + "} {2:>8} {3:^54}").format(
        '%', 'FLUX', 'RXN ID', 'REACTION'))

    for row in consuming.iterrows():
        print_((u"{0.percent:6.1%} {0.flux:>" + str(flux_len) +
                "} {0.id:>8} {0.reaction:>54}").format(row[1]))


def model_summary(model, threshold=1E-8, fva=None, digits=2, **solver_args):
    """Print a summary of the input and output fluxes of the model.

    threshold: float
        tolerance for determining if a flux is zero (not printed)

    fva: int or None
        Whether or not to calculate and report flux variability in the
        output summary

    digits: int
        number of digits after the decimal place to print

    """
    obj_fluxes = pd.Series({'{:<15}'.format(r.id): '{:.3f}'.format(r.x)
                            for r in iterkeys(model.objective)})

    if not fva:

        out_rxns = model.reactions.query(
            lambda rxn: rxn.x > threshold, None
        ).query(lambda x: x, 'boundary')

        in_rxns = model.reactions.query(
            lambda rxn: rxn.x < -threshold, None
        ).query(lambda x: x, 'boundary')

        out_fluxes = pd.Series({r.reactants[0]: r.x for r in out_rxns})
        in_fluxes = pd.Series({r.reactants[0]: r.x for r in in_rxns})

        # sort and round
        out_fluxes.sort_values(ascending=False, inplace=True)
        out_fluxes = out_fluxes.round(digits)
        in_fluxes.sort_values(inplace=True)
        in_fluxes = in_fluxes.round(digits)

        table = pd.np.array(
            [((a if a else ''), (b if b else ''), (c if c else ''))
             for a, b, c in zip_longest(
                ['IN FLUXES'] + in_fluxes.to_string().split('\n'),
                ['OUT FLUXES'] + out_fluxes.to_string().split('\n'),
                ['OBJECTIVES'] + obj_fluxes.to_string().split('\n'))])

    else:
        boundary_reactions = model.reactions.query(lambda x: x, 'boundary')

        fva_results = pd.DataFrame(
            flux_variability_analysis(model, reaction_list=boundary_reactions,
                                      fraction_of_optimum=fva,
                                      **solver_args)).T

        half_span = (fva_results.maximum - fva_results.minimum) / 2
        median = fva_results.minimum + half_span
        rxn_data = pd.concat([median, half_span], 1)
        rxn_data.columns = ['x', 'err']

        for r in rxn_data.index:
            rxn_data.loc[r, 'met'] = model.reactions.get_by_id(r).reactants[0]

        rxn_data.set_index('met', drop=True, inplace=True)

        out_fluxes = rxn_data[rxn_data.x > threshold]
        in_fluxes = rxn_data[rxn_data.x < -threshold]

        out_fluxes = out_fluxes.sort_values(by='x', ascending=False)
        out_fluxes = out_fluxes.round(digits)
        in_fluxes = in_fluxes.sort_values(by='x')
        in_fluxes = in_fluxes.round(digits)

        in_fluxes_s = in_fluxes.apply(
            lambda x: u'{0:0.2f} \u00B1 {1:0.2f}'.format(x.x, x.err),
            axis=1)
        out_fluxes_s = out_fluxes.apply(
            lambda x: u'{0:0.2f} \u00B1 {1:0.2f}'.format(x.x, x.err),
            axis=1)
        out_fluxes_s = out_fluxes.apply(lambda x: text_type(x.x) +
                                        u" \u00B1 " + text_type(x.err), axis=1)

        table = pd.np.array(
            [((a if a else ''), (b if b else ''), (c if c else ''))
             for a, b, c in zip_longest(
                ['IN FLUXES'] + in_fluxes_s.to_string().split('\n'),
                ['OUT FLUXES'] + out_fluxes_s.to_string().split('\n'),
                ['OBJECTIVES'] + obj_fluxes.to_string().split('\n'))])

    print_(u'\n'.join([u"{a:<30}{b:<30}{c:<20}".format(a=a, b=b, c=c) for
                       a, b, c in table]))
