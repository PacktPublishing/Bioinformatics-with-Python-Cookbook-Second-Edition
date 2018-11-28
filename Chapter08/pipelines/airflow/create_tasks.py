""" Task creator
"""

from datetime import datetime, timedelta
import os
from os.path import isfile

from airflow import DAG
from airflow.operators.python_operator import PythonOperator
from airflow.contrib.hooks.ftp_hook import FTPHook
from airflow.contrib.hooks.fs_hook import FSHook
from airflow.contrib.sensors.file_sensor import FileSensor

dag_args = {
    'owner': 'airflow',
    'description': 'Bioinformatics with Python Cookbook pipeline',
    'depends_on_past': False,
    'start_date': datetime(2016, 1, 18),
    'email': ['your@email.here'],
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=1)
}
dag = DAG('bioinf', default_args=dag_args, schedule_interval=None)

ftp_directory = '/hapmap/genotypes/hapmap3/plink_format/draft_2/'
ftp_files = {
    'hapmap3_r2_b36_fwd.consensus.qc.poly.map.bz2': 'hapmap.map.bz2',
    'hapmap3_r2_b36_fwd.consensus.qc.poly.ped.bz2': 'hapmap.ped.bz2',
    'relationships_w_pops_121708.txt': 'relationships.txt',
}


def download_files(ds, **kwargs):
    fs = FSHook('fs_bioinf')
    force = kwargs['params'].get('force', 'false') == 'true'
    with FTPHook('ftp_ncbi') as ftp:
        for ftp_name, local_name in ftp_files.items():
            local_path = fs.get_path() + '/' + local_name
            uncompressed_local_path = local_path[:-4]
            if (isfile(local_path) or isfile(uncompressed_local_path)) and not force:
                continue
            if not isfile(local_name):
                ftp.retrieve_file(ftp_directory + ftp_name, local_path)
    open(fs.get_path() + '/done.txt', 'wb')
    return True

s_fs_sensor = FileSensor(
    task_id='download_sensor',
    fs_conn_id='fs_bioinf',
    filepath='done.txt',
    dag=dag)


def uncompress_files(ds, **kwargs):
    fs = FSHook('fs_bioinf')
    if isfile(fs.get_path() + '/hapmap.ped'):
        return True
    os.system('bzip2 -d {fs_path}/hapmap.map.bz2'.format(fs_path=fs.get_path()))
    os.system('bzip2 -d {fs_path}/hapmap.ped.bz2'.format(fs_path=fs.get_path()))
    return True


def subsample_10p(ds, **kwargs):
    fs = FSHook('fs_bioinf')
    os.system('/home/tra/anaconda3/bin/plink --recode --file {fs_path}/hapmap --noweb --out {fs_path}/hapmap10 --thin 0.1 --geno 0.1'.format(fs_path=fs.get_path()))
    return True


def subsample_1p(ds, **kwargs):
    fs = FSHook('fs_bioinf')
    os.system('/home/tra/anaconda3/bin/plink --recode --file {fs_path}/hapmap --noweb --out {fs_path}/hapmap1 --thin 0.01 --geno 0.1'.format(fs_path=fs.get_path()))
    return True


def compute_pca(ds, **kwargs):
    fs = FSHook('fs_bioinf')
    os.system('/home/tra/anaconda3/bin/plink --pca --file {fs_path}/hapmap1 -out {fs_path}/pca'.format(fs_path=fs.get_path()))
    return True


def plot_pca(ds, **kwargs):
    import matplotlib
    matplotlib.use('svg')
    import pandas as pd
    fs = FSHook('fs_bioinf')
    pca_df = pd.read_csv(fs.get_path() + '/pca.eigenvec', sep=' ', header=None)
    ax = pca_df.plot.scatter(x=2, y=3)
    ax.figure.savefig(fs.get_path() + '/pca.png')


t_dowload_files = PythonOperator(
    task_id='download_files',
    provide_context=True,
    python_callable=download_files,
    params={'force': 'Force download (boolean)'},
    dag=dag)

#Explain below
s_fs_sensor.set_upstream(t_dowload_files)

t_uncompress_files = PythonOperator(
    task_id='uncompress_files',
    provide_context=True,
    python_callable=uncompress_files,
    dag=dag)
t_uncompress_files.set_upstream(s_fs_sensor)

t_subsample_10p = PythonOperator(
    task_id='subsample_10p',
    provide_context=True,
    python_callable=subsample_10p,
    dag=dag)
t_subsample_10p.set_upstream(t_uncompress_files)

t_subsample_1p = PythonOperator(
    task_id='subsample_1p',
    provide_context=True,
    python_callable=subsample_1p,
    dag=dag)
t_subsample_1p.set_upstream(t_uncompress_files)

t_compute_pca = PythonOperator(
    task_id='compute_pca',
    provide_context=True,
    python_callable=compute_pca,
    dag=dag)
t_compute_pca.set_upstream(t_subsample_1p)

t_plot_pca = PythonOperator(
    task_id='plot_pca',
    provide_context=True,
    python_callable=plot_pca,
    dag=dag)
t_plot_pca.set_upstream(t_compute_pca)




