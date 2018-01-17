from enum import Enum


class Group(Enum):
    system = 'System'
    protein = 'Protein'
    protein_hydrogen_less = 'Protein-H'
    alpha_carbons = 'C-alpha'
    backbone = 'Backbone'
    main_chain = 'MainChain'
    main_chain_with_beta_carbons = 'MainChain+Cb'

    def __str__(self):
        return self.value
