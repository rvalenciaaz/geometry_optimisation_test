# download_molecules.py

import sys
import logging
from chembl_webresource_client.new_client import new_client
from sqlalchemy import create_engine, Column, String, Integer, Base
from sqlalchemy.orm import sessionmaker
import selfies as sf

# Set up logging
logging.basicConfig(filename='download_molecules.log', level=logging.INFO,
                    format='%(asctime)s %(levelname)s:%(message)s')

# Database setup
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

class Molecule(Base):
    __tablename__ = 'molecules'
    id = Column(Integer, primary_key=True)
    chembl_id = Column(String, unique=True)
    selfies = Column(String)

def download_molecules(max_records, db_url):
    # Initialize database
    engine = create_engine(db_url)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # Initialize ChEMBL molecule client
    molecule = new_client.molecule

    # Define filters
    filter_criteria = {
        'molecule_properties__mw_freebase__lte': 300,
        'molecule_structures__canonical_smiles__isnull': False,
    }

    # Search molecules
    results = molecule.filter(**filter_criteria).only(['molecule_chembl_id', 'molecule_structures'])

    # Limit number of records
    results = results[:max_records]

    # Process and store molecules
    for res in results:
        chembl_id = res['molecule_chembl_id']
        smiles = res['molecule_structures']['canonical_smiles']
        # Convert SMILES to SELFIES
        try:
            selfies_str = sf.encoder(smiles)
            mol = Molecule(chembl_id=chembl_id, selfies=selfies_str)
            session.add(mol)
            logging.info(f"Successfully added molecule {chembl_id}")
        except Exception as e:
            logging.error(f"Failed to process molecule {chembl_id}: {e}")

    session.commit()
    session.close()

def main():
    max_records = int(sys.argv[1])  # e.g., 1000
    db_url = sys.argv[2]            # e.g., 'sqlite:///molecules.db'
    download_molecules(max_records, db_url)

if __name__ == '__main__':
    main()
