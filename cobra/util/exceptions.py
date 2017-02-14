class DefunctError(Exception):
    """Exception for retired functionality

    Parameters
    ----------
    what : string
        The name of the retired object
    alternative : string
        Suggestion for an alternative
    url : string
        A url to alternative resource
    """

    def __init__(self, what, alternative=None, url=None):
        message = "{} has been removed from cobrapy".format(what)
        if alternative is None:
            message += (" without replacement. Raise an issue at "
                        "https://github.com/opencobracobrapy if you miss it.")
        if alternative is not None:
            message += ". Consider using '{}' instead".format(alternative)
        if url is not None:
            message += " [{}]".format(url)
        super(DefunctError, self).__init__(message)
