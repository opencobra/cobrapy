"""Define the singleton meta class."""


class Singleton(type):
    """Implementation of the singleton pattern as a meta class."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        """Override an inheriting class' call."""
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
