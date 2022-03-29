"""Test the behaviour of the ProcessPool class."""


import os
from typing import Iterable, Tuple

import pytest
from pytest_mock import MockerFixture

from cobra.util import ProcessPool


def dummy_initializer(*args: Iterable) -> Tuple:
    """Implement a 'do nothing' function that accepts initialization arguments."""
    return args


def square(num: int) -> int:
    """Return the square of an integer."""
    return num * num


def summation(*args: Iterable[int]) -> int:
    """Return the sum of all integer arguments."""
    return sum(args)


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
@pytest.mark.parametrize(
    "attributes",
    [
        {},
        {"processes": 2},
        {"initializer": dummy_initializer},
        {"initializer": dummy_initializer, "initargs": (1, "2", [3], {"a": 4})},
        {"maxtasksperchild": 1},
    ],
)
def test_init(attributes: dict) -> None:
    """Test that a process pool can be initialized with each of its arguments."""
    ProcessPool(**attributes)


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_close(mocker: MockerFixture) -> None:
    """Test that the composed pool is closed as well."""
    pool = ProcessPool(processes=3)
    mock = mocker.patch.object(pool, "_pool", autospec=True)
    pool.close()
    mock.close.assert_called_once()


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_with_context(mocker: MockerFixture) -> None:
    """Test that the composed pool's context is managed as well."""
    pool = ProcessPool(processes=3)
    mock = mocker.patch.object(pool, "_pool", autospec=True)
    with pool:
        pass
    mock.__enter__.assert_called_once()
    mock.__exit__.assert_called_once()


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_apply() -> None:
    """Test that a function can be applied."""
    with ProcessPool(processes=3) as pool:
        assert pool.apply(square, (3,)) == 9


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_apply_async() -> None:
    """Test that a function can be applied asynchronously."""
    with ProcessPool(processes=3) as pool:
        assert pool.apply_async(square, (3,)).get() == 9


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_map() -> None:
    """Test that a function can be mapped over an iterable of values."""
    with ProcessPool(processes=3) as pool:
        assert sum(pool.map(square, [2] * 6)) == 24


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_map_async() -> None:
    """Test that a function can be mapped over an iterable of values asynchronously."""
    with ProcessPool(processes=3) as pool:
        assert sum(pool.map_async(square, [2] * 6).get()) == 24


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_imap() -> None:
    """Test that mapped function results can be iterated."""
    with ProcessPool(processes=3) as pool:
        total = 0
        for result in pool.imap(square, [2] * 6):
            total += result
        assert total == 24


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_imap_unordered() -> None:
    """Test that mapped function results can be iterated in any order."""
    with ProcessPool(processes=3) as pool:
        assert sum(pool.imap_unordered(square, [2] * 6)) == 24


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_starmap() -> None:
    """Test that a function can be starmapped over many iterables."""
    with ProcessPool(processes=3) as pool:
        assert sum(pool.starmap(summation, [range(10), range(10), range(10)])) == 135


@pytest.mark.skipif("SKIP_MP" in os.environ, reason="unsafe for parallel execution")
def test_starmap_async() -> None:
    """Test that a function can be starmapped over many iterables asynchronously."""
    with ProcessPool(processes=3) as pool:
        assert (
            sum(pool.starmap_async(summation, [range(10), range(10), range(10)]).get())
            == 135
        )
