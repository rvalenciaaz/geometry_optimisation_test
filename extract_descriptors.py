# extract_descriptors.py

import sys
import logging
import multiprocessing as mp
from sqlalchemy import create_engine, Column, String, Integer, LargeBinary, Float, Base
from sqlalchemy.orm import sessionmaker
from rdkit import Chem
from rdkit.Chem import Descriptors3D

# Set up logging
logging.basicConfig(filename='extract_descriptors.log', level=logging.INFO,
                    format='%(asctime)s %(levelname)s:%(message)s')

# Database setup
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True)
    chembl_id = Column(String, unique=True)
    optimized_geometry = Column(LargeBinary)
    pmi1 = Column(Float)
    pmi2 = Column(Float)
    npr1 = Column(Float)
    npr2 = Column(Float)

def compute_descriptors(mol_data):
    chembl_id, optimized_geometry = mol_data
    # Database session
    engine = create_engine(db_url)
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        xyz_block = optimized_geometry.decode('utf-8')
        mol = Chem.MolFromXYZBlock(xyz_block)
        if mol is None:
            raise ValueError("Invalid molecule from XYZ")
        descriptors = {
            'pmi1': Descriptors3D.PMI1(mol),
            'pmi2': Descriptors3D.PMI2(mol),
            'npr1': Descriptors3D.NPR1(mol),
            'npr2': Descriptors3D.NPR2(mol),
        }
        # Update database
        session.query(Molecule).filter(Molecule.chembl_id == chembl_id).update(descriptors)
        session.commit()
        logging.info(f"Descriptors computed for molecule {chembl_id}")
    except Exception as e:
        logging.error(f"Failed to compute descriptors for molecule {chembl_id}: {e}")
    finally:
        session.close()

def main(db_url, num_processes):
    # Initialize database
    engine = create_engine(db_url)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Fetch molecules without descriptors
    molecules = session.query(Molecule).filter(Molecule.optimized_geometry != None, Molecule.pmi1 == None).all()
    mol_data_list = [(mol.chembl_id, mol.optimized_geometry) for mol in molecules]
    session.close()

    # Set up multiprocessing
    with mp.Pool(processes=num_processes) as pool:
        pool.map(compute_descriptors, mol_data_list)

if __name__ == '__main__':
    db_url = sys.argv[1]         # e.g., 'sqlite:///molecules.db'
    num_processes = int(sys.argv[2])  # Number of parallel processes
    main(db_url, num_processes)
