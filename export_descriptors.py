# export_descriptors.py

import sys
from sqlalchemy import create_engine
import pandas as pd

def export_descriptors(db_url, output_file):
    engine = create_engine(db_url)
    query = "SELECT chembl_id, pmi1, pmi2, npr1, npr2 FROM molecules WHERE pmi1 IS NOT NULL"
    df = pd.read_sql(query, engine)
    df.to_csv(output_file, index=False)

if __name__ == '__main__':
    db_url = sys.argv[1]         # e.g., 'sqlite:///molecules.db'
    output_file = sys.argv[2]    # e.g., 'descriptors.csv'
    export_descriptors(db_url, output_file)
