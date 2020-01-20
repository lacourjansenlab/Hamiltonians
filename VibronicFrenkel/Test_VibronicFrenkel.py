import unittest
import numpy as np

class TestFranckCondon(unittest.TestCase):
    def setUp(self):
        from VibronicFrenkel import Autovar
        parameters = {
            "N": 21,
            "E": 0,
            "J": 7.9,
            "MaxVib": 1,
            "wvib": 1110.12,
            "S": 0.774,
            "sigma": 0,
            "units": "wavenumbers"
        }
        self.p = Autovar(parameters)

    def test_franck_condon(self):
        from VibronicFrenkel import gen_fc
        # Derivation can be found in
        # Modern Optical Spectroscopy with Examples from Biophysics and Biochemistry
        # by: W.W. Parson
        # page: 167

        # Check 0-0
        fc = gen_fc(0,0,self.p.S)
        self.assertEqual(fc,np.exp(-self.p.S/2))

        # Check 0-1
        fc = gen_fc(0,1,self.p.S)
        zero_one = np.exp(-self.p.S/2)*np.sqrt(self.p.S)
        self.assertTrue(np.allclose(fc,zero_one,rtol=1e-5,atol=1e-8))

        # Check 3-1
        fc = gen_fc(3,1,self.p.S)
        two_zero = np.exp(-self.p.S/2)*self.p.S/np.sqrt(2)
        three_zero = np.exp(-self.p.S/2)*(-1)**3*self.p.S**(3/2)/np.sqrt(6)
        three_one = (np.sqrt(3)*two_zero+np.sqrt(self.p.S)*three_zero)
        self.assertTrue(np.allclose(fc,three_one,rtol=1e-5,atol=1e-8))

# S (Huang-Rhys factor) = 0 testcase
class TestHamiltonian_ZeroHR(unittest.TestCase):
    def setUp(self):
        from VibronicFrenkel import Autovar
        parameters = {
            "N": 5,
            "E": 0,
            "J": 786.9,
            "MaxVib": 1,
            "wvib": 1110.12,
            "S": 0,
            "sigma": 0,
            "units": "wavenumbers"
        }
        self.p = Autovar(parameters)

        # Generate Franck-Condon factors between ground state and excited state
        from VibronicFrenkel import gen_fc_table
        fc_table = gen_fc_table(self.p.S,self.p.MaxVib)

        # Index one-particle states
        from VibronicFrenkel import index_ops
        ops = index_ops(self.p)

        from VibronicFrenkel import gen_hamiltonian
        self.H = gen_hamiltonian(self.p,ops,fc_table)

    def test_size(self):
        ops_size = self.p.N*(self.p.MaxVib+1)
        self.assertEqual(np.shape(self.H),(ops_size,ops_size))

    def test_eigenvalues(self):
        # Eigenvalues should be same as closed boundary condition problems.
        # But with the vibrational states completely decoupled
        Nspace = np.linspace(1,self.p.N,self.p.N)

        # Theoretical eigenvalues
        Eig_t_ground = self.p.E+2*self.p.J*np.cos(np.pi*Nspace/(self.p.N+1))
        Eig_t_vib = np.zeros((self.p.N,)) + self.p.E + self.p.wvib
        Eig_t = np.concatenate((Eig_t_ground,Eig_t_vib))
        Eig_t = np.sort(Eig_t)

        # Numerical eigenvalues
        [Eig,Vec] = np.linalg.eig(self.H)
        Eig = np.sort(Eig)

        # Compare
        self.assertTrue(np.allclose(Eig,Eig_t))

if __name__ == '__main__':
    unittest.main()
