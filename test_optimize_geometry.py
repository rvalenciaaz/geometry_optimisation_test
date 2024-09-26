# test_optimize_geometry.py

import unittest
from optimize_geometry import optimize_molecule

class TestOptimizeGeometry(unittest.TestCase):
    def test_optimize_valid_molecule(self):
        # Example SELFIES string for water molecule
        selfies_str = '[O][H][H]'
        mol_data = ('test_molecule', selfies_str)
        # Set up a mock database URL
        global db_url
        db_url = 'sqlite:///:memory:'
        optimize_molecule(mol_data)
        # Assertions can be added here to check database entries

if __name__ == '__main__':
    unittest.main()
