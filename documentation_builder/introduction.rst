Getting Started with cobrapy
============================

To begin with, cobrapy comes with two bundled models for Salmonella and E. coli.
To load a test model, type

>>> import cobra.test
>>> model = cobra.test.create_test_model(cobra.test.salmonella_pickle)

The reactions, metabolites, and genes attributes of the cobrapy model are 
are a special type of list called a :class:`~cobra.core.DictList.DictList`, 
and each one is made up of :class:`~cobra.core.Reaction.Reaction`, 
:class:`~cobra.core.Metabolite.Metabolite` and 
:class:`~cobra.core.Gene.Gene` objects respectively.

>>> print len(model.reactions)
2546
>>> print len(model.metabolites)
1802
>>> print len(model.genes)
1271

Just like a regular list, objects in the DictList can be retrived by index. 
For example, to get the 30th reaction in the model:

>>> print model.reactions[29]
2AGPA180tipp

Addictionally, items can be retrived by their id using the
:func:`~cobra.core.DictList.DictList.get_by_id` function. For example, to get
the cytosolic atp metabolite (the id is "atp_c"), we can do the following:

>>> atp_c = model.metabolites.get_by_id("atp_c")
>>> print atp_c.formula
C10H12N5O13P3

In the above example, atp_c is a :class:`~cobra.core.Metabolite.Metabolite`
object, which has a formula attribute.

As an added bonus, users with an interactive shell such as IPython will be able
to tab-complete to list elements inside a list. While this is not recommended
behavior for most code because of the possibility for characters like "-" inside
ids, this is very useful while in an interactive prompt:

>>> model.reactions.EX_glc__D_e.lower_bound
0.0
