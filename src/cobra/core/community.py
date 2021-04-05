"""Define Community Model class."""

from copy import deepcopy

from cobra.core.dictlist import DictList
from cobra.core.model import Model


class Community(Model):
    """Class representation of cobra community model.

    Parameters
    ----------
    models : list
        List of model objects.
    id_ : string
        Model identifier.
    name : string
        Human readable name for the model.
    """

    def __init__(self, models, id_="community_model", name=None):
        """Init community model."""
        Model.__init__(self, id_, name)

        self.organisms = DictList()

        exchange_reactions = DictList()

        for model in models:
            m = deepcopy(model)

            exchange_reactions.union(deepcopy(m.exchanges))

            for metabolite in m.metabolites:
                metabolite.id = f"{metabolite.id}_{model.id}"

            for gene in m.genes:
                gene.id = f"{gene.id}_{model.id}"

            for reaction in m.reactions:
                if reaction.id.startswith("EX_"):
                    exchange = exchange_reactions.get_by_id(reaction.id)
                    reaction -= exchange
                reaction.id = f"{reaction.id}_{model.id}"

                # translate gpr
                gpr = reaction.gene_reaction_rule
                gpr = gpr.replace("(", "( ")
                gpr = gpr.replace(")", " )")
                reaction.gene_reaction_rule = " ".join(
                    [
                        f"{token}_{model.id}"
                        if token not in ["and", "or", "(", ")"]
                        else token
                        for token in gpr.split()
                    ]
                )

            for group in m.groups:
                group.id = f"{group.id}_{model.id}"

            self.add_metabolites(m.metabolites)
            self.add_reactions(m.reactions)
            self.add_groups(m.groups)

            if self.objective.expression:
                self.objective = self.objective.expression + m.objective.expression
            else:
                self.objective = m.objective

            self.organisms.append(m)

        self.add_reactions(exchange_reactions)
