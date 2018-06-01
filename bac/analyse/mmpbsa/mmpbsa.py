import os
import subprocess


class MMPBSA:

    def __init__(self):
        self._salt_concentration = 0.0

    @property
    def salt_concentration(self):
        return self._salt_concentration

    @salt_concentration.setter
    def salt_concentration(self, value):
        self._salt_concentration = value

    def get_script(self):

        return """&general
startframe=1, endframe=400, interval=100, keep_files=0,
verbose=2, keep_files=2,
strip_mask=':WAT,Cl*,CIO,Cs+,IB,K*,Li+,MG*,Na+,Rb+,CS,RB,NA,F,CL,ZN'
/
&gb
igb=5, saltcon={},
/
&pb
istrng=0.0, fillratio=4.0, inp=1, radiopt=0
/
        """.format(self.salt_concentration)

    def get_free_energy(self, trajectory, solvated_complex, dry_complex, receptor, ligand, output):

            if 'AMBERHOME' not in os.environ.keys():
                raise EnvironmentError("Set $AMBERHOME for this to work!")

            with open('mmpbsa.in', 'w') as f:
                f.write(self.get_script())

            executable = f"MMPBSA.py -O -i mmpbsa.in -o {output} -sp {solvated_complex} " \
                         f"-cp {dry_complex} -rp {receptor} -lp {ligand} -y {trajectory}"

            output = subprocess.check_output(executable.split()).decode()

            print(output)
