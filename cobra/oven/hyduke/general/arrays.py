from numpy import ndarray
class ResultsArray(ndarray):
    """A primitive wrapper to allow accessing numpy.ndarrays via
    named rows and columns.  The ResultsArray.row_names and
    column_names must be assigned after the object is created.

    The names will not carry over for any operations.

    TODO: Finish the implementation

    """
    def __init__(self, shape, row_names=None, column_names=None):
        ndarray.__init__(shape)
        if row_names:
            self.row_names = row_names
        else:
            self.row_names = range(shape[0])
        if column_names:
            self.column_names = column_names
        else:
            column_names = range(shape[1])
    def get(self, row_name=None, column_name=None):
        if row_name:
            the_row = self.row_names.index(row_name)
        if column_name:
            the_column = self.column_names.index(column_name)
        if row_name and column_name:
            return self[the_row, the_column]
        if not row_name:
            return self[:, the_column]
        if not column_name:
                return self[the_row, :]
