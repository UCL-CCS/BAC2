import os
import yaml
from enum import Enum
import warnings

import parmed as pmd

from bac.utils.decorators import pathlike, advanced_property

#Defaults
DEFAULT_INFO_PATH = os.path.dirname(os.path.realpath(__file__))

class FFType(Enum):
    amber = 'amber'

class SystemBuilder():

    def __init__(self, ff_type=FFType('amber'), system_type='protein_small_ligand', ff=['ff14sb', 'tip3p'], sanitize=False, **kwargs):

        self.ff_type = ff_type
        self.system_type = system_type
        self.ff_elements = ff

        self.system_path = kwargs.get('system_path', DEFAULT_INFO_PATH)
        self.ff_add = kwargs.get('ff_add', [])

        self.ingredients = kwargs.get('ingredients', {})

        description_path = os.path.join(self.system_path, 'system_descriptions', system_type + '.yml')

        if not os.path.isfile(description_path):
            raise IOError(f"System description file does not exist: {description_path}")

        stream = open(description_path, 'r')
        self.description = yaml.load(stream)

        expected_ingredients = self.description.keys()

        if self.ingredients:
            actual_ingredients = self.ingredients.keys()
            if actual_ingredients != expected_ingredients:
                warnings.warn(f"Expected ingredients {expected_ingredients} but provided {actual_ingredients}")
        else:
            for ingredient_type in self.description.keys():
                self.ingredients[ingredient_type] = []


    @advanced_property(type=FFType, default=FFType.amber)
    def ff_type(self): pass

    @property
    def system_type(self): pass

    @advanced_property(type=list, default=['ff14sb', 'tip3p'])
    def ff_elements(self): pass

    @advanced_property(type=dict, default={})
    def ingredients(self): pass

    @advanced_property(type=dict, default={})
    def description(self): pass

    def needed_ingredients(self):

        needed = []

        for type, files in self.ingredients.items():
            if not files:
                needed.append(type)

        return needed

    def needed_ff_moltypes(self):

        for

        pass

    def check_ingredients(self):

        pass
