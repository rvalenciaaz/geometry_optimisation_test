python download_molecules.py 1000 sqlite:///molecules.db

python optimize_geometry.py sqlite:///molecules.db 4

python extract_descriptors.py sqlite:///molecules.db 4


# Build the Docker image
docker build -t molecule_pipeline .

# Run the Docker container
docker run -p 8080:8080 molecule_pipeline


python -m unittest test_optimize_geometry.py

python export_descriptors.py sqlite:///molecules.db descriptors.csv


docker build -t molecule_pipeline .
docker run -p 8080:8080 molecule_pipeline


