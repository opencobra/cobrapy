"""Test functions of context.py."""

from typing import TYPE_CHECKING

import pytest

from cobra.util import HistoryManager, get_context, resettable


if TYPE_CHECKING:
    from cobra import Model


def test_history_manager() -> None:
    """Test initialization and resetting of HistoryManager."""
    # initialize HistoryManager
    history_manager = HistoryManager()
    # add non-functioning operation
    history_manager(lambda: 1)
    # reset operations
    history_manager.reset()


def test_get_context(model: "Model") -> None:
    """Test if context retrieval is working."""
    with model as model:
        # reverse optimization direcion
        model.objective_direction = "min"
        context = get_context(model)
        # should have history
        assert context._history[-1].args[-1] == "max"

    # there shouldn't be any history
    with pytest.raises(IndexError):
        context._history[-1]


def test_resettable() -> None:
    """Test if resettable decorator is functional."""
    # decorate a dummy function
    @resettable
    def change_my_name(old_name, new_name):
        """Change old name to new name."""
        if old_name != new_name:
            old_name = new_name

    # call the dummy function
    change_my_name("hmm", "hmmm")
