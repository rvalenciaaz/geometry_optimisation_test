# optimize_geometry.py

import os
import sys
import logging
import multiprocessing as mp
from sqlalchemy import create_engine, Column, String, Integer, LargeBinary, Base
from sqlalchemy.orm import sessionmaker
from rdkit import Chem
from rdkit.Chem import AllChem
import psi4
import selfies as sf

# Set up logging
logging.basicConfig(filename='optimize_geometry.log', level=logging.INFO,
                    format='%(asctime)s %(levelname)s:%(message)s')

# Database setup
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True)
    chembl_id = Column(String, unique=True)
    selfies = Column(String)
    optimized_geometry = Column(LargeBinary)

def optimize_molecule(mol_data):
    chembl_id, selfies_str = mol_data
    # Database session
    engine = create_engine(db_url)
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        # Decode SELFIES to SMILES
        smiles = sf.decoder(selfies_str)
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        # Generate XYZ coordinates
        xyz = Chem.MolToXYZBlock(mol)
        # Set Psi4 options
        psi4.core.be_quiet()
        psi4.set_options({
            'reference': 'rhf',
            'basis': 'sto-3g',
            'maxiter': 200,
        })
        # Define geometry
        geometry = psi4.geometry(xyz)
        # Optimize geometry
        psi4.optimize('scf', molecule=geometry)
        # Save optimized geometry
        optimized_xyz = geometry.save_string_xyz()
        # Update database
        session.query(Molecule).filter(Molecule.chembl_id == chembl_id).update(
            {"optimized_geometry": optimized_xyz.encode('utf-8')})
        session.commit()
        logging.info(f"Successfully optimized molecule {chembl_id}")
    except Exception as e:
        logging.error(f"Failed to optimize molecule {chembl_id}: {e}")
    finally:
        session.close()

def main(db_url, num_processes):
    # Initialize database
    engine = create_engine(db_url)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Fetch molecules without optimized geometry
    molecules = session.query(Molecule).filter(Molecule.optimized_geometry == None).all()
    mol_data_list = [(mol.chembl_id, mol.selfies) for mol in molecules]
    session.close()

    # Set up multiprocessing
    with mp.Pool(processes=num_processes) as pool:
        pool.map(optimize_molecule, mol_data_list)

if __name__ == '__main__':
    db_url = sys.argv[1]         # e.g., 'sqlite:///molecules.db'
    num_processes = int(sys.argv[2])  # Number of parallel processes
    main(db_url, num_processes)