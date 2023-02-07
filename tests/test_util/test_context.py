"""Test functions of context.py."""

from typing import TYPE_CHECKING

from cobra.util import HistoryManager, get_context, resettable


if TYPE_CHECKING:
    from cobra import Model


def test_history_manager() -> None:
    """Test initialization and resetting of HistoryManager."""
    # initialize HistoryManager
    history_manager = HistoryManager()
    # add non-functioning operation
    history_manager(lambda: 1)
    # check size of the stack
    assert history_manager.size() == 1
    # reset operations
    history_manager.reset()


def test_get_context(model: "Model") -> None:
    """Test if context retrieval is working."""
    with model as model:
        # reverse optimization direcion
        model.objective_direction = "min"
        # retrieve context
        context = get_context(model)
        # check size of the context
        if context:
            assert context.size() == 1

    # there shouldn't be any history
    if context:
        assert context.size() == 0


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
