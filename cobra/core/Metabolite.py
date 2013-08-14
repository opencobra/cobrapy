#cobra.core.Metabolite.py
#######################
#BEGIN Class Metabolite
#
#Is it better to restrict a Metabolite to a single model or
#should we allow a Metabolite to be associated with multiple models?
#
from .Species import Species
class Metabolite(Species):
    """Metabolite is a class for holding information regarding
    a metabolite in a cobra.Reaction object.

        
    """
    ## __slots__ = ['id', 'formula', 'name', 'compartment',
    ##              'charge']
    ##def __setstate__(self, the_dict):
    ##    from cobra.core.Metabolite import Metabolite
    ##    Species.__setstate__(self, the_dict)
    ##    [self.__setattr__(k, v) for k, v in the_dict]

    def __init__(self, id=None, formula=None,
                 name=None, compartment=None):
        """
        id: A string.

        formula: cobra.Formula or String  of a chemical formula.  Defaults to None
        to save time in pickling and such.
        
        name: String.  A human readable name.
        
        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object
        
        """
        Species.__init__(self, id, formula, name, compartment)
        self._constraint_sense = 'E'
        self._bound = 0.



    def remove_from_model(self, model=None, method='subtractive'):
        """Removes the association

        model: None or :class:`~cobra.core.Model` object.
            Remove the reaction from this model.

        method: 'subtractive' or 'destructive'.  If 'subtractive' then
        the metabolite is removed from all associated reactions.  If
        'destructive' then all associated reactions are removed from the Model.

    
        
        """
        # why is model being taken in as a parameter? This plays
        #back to the question of allowing a Metabolite to be associated
        #with multiple Models
        if model != self._model and model is not None:
            raise Exception('%s not in %s ergo it cannot be removed. (%s)'%(self,
                                                                  model,
                                                                  self._model))
                                                            
        self._model.metabolites.remove(self)
        self._model = None
        if method.lower() == 'subtractive':
            for the_reaction in list(self._reaction):
                the_coefficient = the_reaction._metabolites[self]
                the_reaction.subtract_metabolites({self: the_coefficient})
        elif method.lower() == 'destructive':
            [x.remove_from_model() for x in self._reaction()]
        else:
            raise Exception("method for remove_from_model must be 'subtractive' " +\
                            "or 'destructive'.  You entered: '%s'"%method)

#
#END Class Metabolite
########################
