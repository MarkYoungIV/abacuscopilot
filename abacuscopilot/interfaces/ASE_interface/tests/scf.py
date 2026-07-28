import tempfile
import unittest
from pathlib import Path

here = Path(__file__).parent
from abacuslite import Abacus, AbacusProfile
from abacuslite.io.generalio import load_orbital, load_pseudo
from ase.build import bulk


class TestSCF(unittest.TestCase):

    def test(self):
        pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

        silicon = bulk('Si', 'diamond', a=5.43)
        aprof = AbacusProfile(
            command='mpirun -np 2 abacus',
            pseudo_dir=pporb,
            orbital_dir=pporb,
            omp_num_threads=1,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            abacus = Abacus(
                profile=aprof,
                directory=tmpdir,
                pseudopotentials=load_pseudo(pporb),
                basissets=load_orbital(pporb, efficiency=True),
                inp={
                    'basis_type': 'lcao',
                    'gamma_only': True,
                    'scf_thr': 1e-3, # fast for test, wrong for production
                }
            )

            silicon.calc = abacus
            print('Silicon :', silicon.get_potential_energy())

if __name__ == '__main__':
    unittest.main()
