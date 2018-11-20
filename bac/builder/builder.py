import os
import glob
import yaml
from enum import Enum
import warnings

import parmed as pmd

from supproperty import supproperty

#Defaults
DEFAULT_INFO_PATH = os.path.dirname(os.path.realpath(__file__))

class FFType(Enum):
    amber = 'amber'

class SystemBuilder():

    def __init__(self, ff_type=FFType('amber'), system_type='protein_small_ligand', ff=['ff14sb', 'tip3p'], sanitize=False, **kwargs):

        self.ff_type = ff_type
        self.system_type = system_type

        self.system_path = kwargs.get('system_path', DEFAULT_INFO_PATH)

        ff_dir = os.path.join(self.system_path, self.ff_type.value, 'ff')
        ff_search_txt = os.path.join(ff_dir, '*.yml')
        ff_files = glob.glob(ff_search_txt)

        valid_ff = []
        for ff_file in ff_files:
            ff_name = os.path.splitext(os.path.basename(ff_file))[0]
            valid_ff.append(ff_name)

        self.valid_ff = valid_ff

        ff_contents = []

        for ff_element in ff:
            if ff_element in valid_ff:
                ff_location = os.path.join(ff_dir, ff_element +'.yml')
                stream = open(ff_location)
                ff_desc = yaml.load(stream)
                ff_contents += ff_desc['reside_names']
                self.ff_elements.append(ff_element)
            else:
                warnings.warn("{ff_element} not a valid forcefield suggestion")

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

    @supproperty(type=FFType, default=FFType.amber)
    def ff_type(self): pass

    @property
    def system_type(self): pass

    @supproperty(type=list, default=[])
    def ff_elements(self): pass

    @supproperty(type=dict, default={})
    def ingredients(self): pass

    @supproperty(type=dict, default={})
    def description(self): pass

    def needed_ingredients(self):

        needed = []

        for type, files in self.ingredients.items():
            if not files:
                needed.append(type)

        return needed

    def needed_ff_moltypes(self):

        pass

    def check_ingredients(self):

        pass
