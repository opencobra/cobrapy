import libsbml
import collections
from cobra.core.object import Object

QUALIFIER_TYPES = (
    "is", "hasPart", "isPartOf", "isVersionOf", "hasVersion",
    "isHomologTo", "isDescribedBy", "isEncodedBy", "encodes",
    "occursIn", "hasProperty", "isPropertyOf", "hasTaxon",
    "unknown", "bqm_is", "bqm_isDescribedBy", "bqm_isDerivedFrom",
    "bqm_isInstanceOf", "bqm_hasInstance", "bqm_unknown",
)

class Annotation:
    """Annotation is a class used to hold the annotation
    information of each component of SBML model derived from SBase

    """

    def __init__(self):
        self._cvTerms = []
        self.history = History()
        self.listofKeyValue = []




    class MyDict(dict):
        def __setitem__(self, key, value):
            if key == 'resource':
                if isinstance(value, list):
                    for item in value:
                        if not isinstance(item, list) or len(item) != 2:
                            raise TypeError("Each resource must be of the "
                                            "form ['resource', 'provider']")

                    super().__setitem__(key, value)
                else:
                    raise TypeError("All the resources must be present inside a list")
            elif key == "cvTerm":
                if not isinstance(value, CVTerm):
                    raise TypeError("value must be of the form CVTerm")
                else:
                    super().__setitem__(key, value)
            else:
                raise ValueError("Can't accept keys other than 'resource' and 'cvTerm'")


    class CVTerm:

        def __init__(self):
            self._qualifier = ""
            self._groupResource = []

        @property
        def qualifier(self):
            return getattr(self, "_qualifier", None)

        @qualifier.setter
        def qualifier(self, value):
            if value == self.qualifier:
                pass
            elif not isinstance(value, string_types):
                raise TypeError("qualifier must be a string")
            elif value not in QUALIFIER_TYPES:
                raise ValueError("%s is not a valid qualifier", qualifier)
            else:
                self._qualifier = value

        @property
        def groupResource(self):
            return getattr(self, "_groupResource")


    class History:

        def __init__(self):
            self.creator = {
                "first_name" : "",
                "last_name" : "",
                "email" : "",
                "organization_name" : ""
            }
            self._created = libsbml.Date()
            self._modified = []

        @property
        def created(self):
            return self._created.getDateAsString()

        @created.setter
        def created(self, value):
            if not isinstance(value, str):
                raise TypeError("date passed must be a string")
            else:
                result = self._created.setDateAsString(value)
                if result != 0:
                    raise ValueError(libsbml.OperationReturnValue_toString(result))

        @property
        def modified(self):
            return getattr(self, "_modified", [])

        @modified.setter
        def modified(self, value):
            if not isinstance(value, list):
                raise TypeError("value passed must be a list")
            else:
                date = libsbml.Date()
                modified_list = []
                for item in value:
                    if isinstance(item, str):
                        result = date.setDateAsString(item)
                        if result != 0:
                            raise ValueError(libsbml.OperationReturnValue_toString(result))
                        modified_list.append(date.getDateAsString())
                    else:
                        raise TypeError("%s must be of type string")
                self._modified = modified_list

        def add_modified_date(self, date):
            if not isinstance(date, str):
                raise TypeError("date passed must be a string")
            lib_date = libsbml.Date()
            result = lib_date.setDateAsString(date)
            if result != 0:
                raise ValueError(libsbml.OperationReturnValue_toString(result))
            self._modified.append(lib_date.getDateAsString())



    class KeyValuePairs(Object):

        def __init__(self, id=None, name=None):
            Object.__init__(self, id, name)
            self._key = ''
            self._value = ''
            self.uri = ''

        @property
        def key(self):
            return getattr(self, "_key", None)

        @key.setter
        def key(self, inKey):
            if not isinstance(inKey, string_types):
                raise TypeError("the key must be of type string")
            else:
                self._key = inKey

        @property
        def value(self):
            return getattr(self, "_key", None)

        @value.setter
        def value(self, inValue):
            if not isinstance(inValue, string_types):
                raise TypeError("the value must be of type string")
            else:
                self._value = inValue
