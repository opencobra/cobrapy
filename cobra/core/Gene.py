#cobra.Gene.py
#######################
#BEGIN Class Gene
#
from Metabolite import Metabolite
class Gene(Metabolite):
    """A Gene is a special class of metabolite.
    

    TODO: Make design decisions about TUs and such
    """
    def __init__(self, id, formula=None,
                 name=None, compartment=None, strand='+',
                 locus_start=0, locus_end=0, functional=True):
        """
        id: A string.

        formula: cobra.Formula or String  of a chemical formula.  Defaults to None
        to save time in pickling and such.
        
        name: String.  A human readable name.
        
        compartment: None or a dictionary indicating the cellular location
        of the metabolite.  Used when in a cobra.Reaction or Model
        object

        strand: '+' or '-' to indicate forward or reverse strand for DNA.

        locus_start: Int.  The index of the starting base for the gene.

        locus_end: Int. The index of the last base for the gene.

        functional: Boolean.  Indicate whether the gene is functional.  If it
        is not functional then it cannot be used in an enzyme complex nor
        can its products be used.
        
        """
        Metabolite.__init__(self, id, formula=formula,
                            name=name, compartment=compartment)
        self.locus_start = locus_start
        self.locus_end = locus_end
        self.strand = strand
        self.functional = functional
#
#END Class Gene
########################
