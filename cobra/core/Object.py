from ..external.six import iteritems


class Object(object):
    """Defines common behavior of object in cobra.core"""

    def __init__(self, id=None, mnx_id=None):
        """
        id: None or a string

        mnx_id: None or a String of the MetaNetX.org ID for the Object
        """
        self.id = id
        self.mnx_id = mnx_id
        # The following two fields will eventually
        # be objects that enforce basic rules about
        # formatting notes and annotation
        self.notes = {}
        self.annotation = {}

    def __getstate__(self):
        """To prevent excessive replication during deepcopy."""
        state = self.__dict__.copy()
        if '_model' in state:
            state['_model'] = None
        return state

    def guided_copy(self):
        """.. deprecated :: 0.3 use copy direclty"""
        the_copy = self.__class__(self.id)
        for k, v in iteritems(self.__dict__):
            # Don't try to set properties
            if not isinstance(getattr(type(self), k, None), property):
                setattr(the_copy, k, v)
        return(the_copy)

    def _copy_parent_attributes(self, gene_object):
        """Helper function for shallow copying attributes from a parent object
        into a new child object.

        """
        for k, v in iteritems(gene_object.__dict__):
            setattr(self, k, v)

    def startswith(self, x):
        return self.id.startswith(x)

    def endswith(self, x):
        return self.id.endswith(x)

    def __contains__(self, x):
        return self.id.__contains__(x)

    def __getitem__(self, index):
        return self.id[index]

    def __getslice__(self, i, j):
        return self.id[i:j]

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

    def __str__(self):
        return str(self.id)
