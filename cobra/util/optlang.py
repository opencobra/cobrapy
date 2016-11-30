def split_bounds(lower_bound, upper_bound):
    """Split a given (lower_bound, upper_bound) interval into a negative
    component and a positive component. Negative components are negated
    (returns positive ranges) and flipped for usage with forward and reverse
    reactions bounds

    """

    # assert lower_bound <= upper_bound, "lower bound is greater than upper"

    bounds_list = [0, 0, lower_bound, upper_bound]
    bounds_list.sort()

    return -bounds_list[1], -bounds_list[0], bounds_list[2], bounds_list[3]


def update_forward_and_reverse_bounds(reaction):
    """For the given reaction, update the bounds in the forward and
    reverse variable bounds.

    Parameters
    ----------
    reaction : cobra.Reaction object

    """

    r_lb, r_ub, f_lb, f_ub = split_bounds(*reaction.bounds)

    try:
        reaction.reverse_variable.lb = r_lb
        reaction.reverse_variable.ub = r_ub

        reaction.forward_variable.lb = f_lb
        reaction.forward_variable.ub = f_ub

    except AttributeError:
        pass
