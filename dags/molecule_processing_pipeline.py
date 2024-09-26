# molecule_processing_pipeline.py

from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.utils.task_group import TaskGroup
from datetime import datetime, timedelta

default_args = {
    'owner': 'airflow',
    'start_date': datetime(2023, 10, 1),
    'retries': 2,
    'retry_delay': timedelta(minutes=5),
}

with DAG(
    'molecule_processing_pipeline',
    default_args=default_args,
    schedule_interval=None,
    catchup=False,
) as dag:

    download_task = BashOperator(
        task_id='download_molecules',
        bash_command='python /app/download_molecules.py 1000 sqlite:///app/molecules.db',
    )

    optimize_task = BashOperator(
        task_id='optimize_geometry',
        bash_command='python /app/optimize_geometry.py sqlite:///app/molecules.db 4',
    )

    extract_task = BashOperator(
        task_id='extract_descriptors',
        bash_command='python /app/extract_descriptors.py sqlite:///app/molecules.db 4',
    )

    download_task >> optimize_task >> extract_task
