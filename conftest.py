def pytest_addoption(parser):
    parser.addoption("--run-slow", action="store_true",
                     help="run slow tests")
    parser.addoption("--run-non-deterministic", action="store_true",
                     help="run tests that not always yield the same results")
