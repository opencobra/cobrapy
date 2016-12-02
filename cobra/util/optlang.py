def separate_forward_and_reverse_bounds(lower_bound, upper_bound):
    """Split a given (lower_bound, upper_bound) interval into a negative
    component and a positive component. Negative components are negated
    (returns positive ranges) and flipped for usage with forward and reverse
    reactions bounds

    """

    assert lower_bound <= upper_bound, "lower bound is greater than upper"

    bounds_list = [0, 0, lower_bound, upper_bound]
    bounds_list.sort()

    return -bounds_list[1], -bounds_list[0], bounds_list[2], bounds_list[3]


def update_forward_and_reverse_bounds(reaction, direction='both'):
    """For the given reaction, update the bounds in the forward and
    reverse variable bounds.

    Parameters
    ----------
    reaction : cobra.Reaction object

    """

    reverse_lb, reverse_ub, forward_lb, forward_ub = \
        separate_forward_and_reverse_bounds(*reaction.bounds)

    try:
        # Clear the original bounds to avoid complaints
        if direction == 'both':
            reaction.forward_variable._ub = None
            reaction.reverse_variable._lb = None
            reaction.reverse_variable._ub = None
            reaction.forward_variable._lb = None

        if direction in {'both', 'upper'}:
            reaction.forward_variable.ub = forward_ub
            reaction.reverse_variable.lb = reverse_lb

        if direction in {'both', 'lower'}:
            reaction.reverse_variable.ub = reverse_ub
            reaction.forward_variable.lb = forward_lb

    except AttributeError:
        pass
