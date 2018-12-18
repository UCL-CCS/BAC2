"""
Create models ready for parameterization.
"""
import os
import tempfile
import subprocess
import glob
import parmed as pmd
from .source_structure import SourceStructure
from .modeller.loop_model import create_loop_model_script, write_ali_file
from bac.utils.file_system import cd


class ScaffoldBuilder:

    def __init__(self, pdb_filename,
                 modelling_method='modeller',
                 work_dir='.'):

        if modelling_method != 'modeller':
            raise NotImplementedError('At present model completion is only '
                                      'implement using "modeller" ')
        self.modelling_method = modelling_method

        self.source = SourceStructure(pdb_filename)
        self.work_dir = work_dir

    def build_scaffold(self):

        source = self.source

        components = []

        for chain, subdivisions in source.decomposition.items():

            for subdivision in subdivisions:

                initial_structure = subdivision.structure()

                if (subdivision.is_complete() and
                        subdivision.residue_type in ['protein']):

                    if subdivision.gaps:

                        tmp_dir = tempfile.mkdtemp(prefix='scaffold_',
                                                   dir=self.work_dir)

                        structure_id = 'tmp'
                        pdb_filename = f'{structure_id}.pdb'
                        pdb_path = os.path.join(tmp_dir, pdb_filename)

                        # TODO: AltLocs?
                        initial_structure.write_pdb(pdb_path, renumber=False)

                        added_residues = subdivision.get_added_residue_positions()

                        script_name = create_loop_model_script(pdb_filename,
                                                               added_residues,
                                                               workdir=tmp_dir)

                        # Write alignment file
                        target_sequence = subdivision.target_sequence()
                        # Need to think about pre/post residues here
                        stucture_sequence = subdivision.sequence()

                        start_resid = initial_structure.residues[0].number
                        chain_id = initial_structure.residues[0].chain
                        end_resid = initial_structure.residues[-1].number

                        ali_filename = 'alignment.ali'
                        alignment_path = os.path.join(tmp_dir, ali_filename)
                        write_ali_file(alignment_path, target_sequence,
                                       stucture_sequence,
                                       structure_id=structure_id,
                                       start_resid=start_resid,
                                       end_resid=end_resid,
                                       chain=chain_id)

                        # Run script
                        with cd(tmp_dir):

                            subprocess.call(['python', script_name])
                            completed = glob.glob(f'{structure_id}_fill*.pdb')[0]
                            struct = pmd.load_file(completed)

                        # Read new model into structure and add to components
                        components.append(struct)

                    else:

                        components.append(initial_structure)

                elif subdivision.residue_type in ['ions', 'water']:

                    components.append(initial_structure)

        return components
