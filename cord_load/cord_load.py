"""
HOW TO USE
1. Unpack the archive and put the script inside like this:
 ├── cord_load.py
 ├── metadata.csv
 └── document_parses
     └── pdf_json
         ├─ 27f774dd62b31a5a4eea.json
         └── ...
2. Check DB credentials at DB_AUTH_STRING (find string below)
3. Launch the script, wait for 'Done' message in the log

Call create_stripped_archive() to strip inital CORD-data
"""

import os
from shutil import copyfile
import time
from uuid import uuid4
from datetime import datetime
import pandas as pd
import psycopg2
import psycopg2.extras as extras

class CordParser:
    """
    Processing CORD files: loading, parsing and inserting to the DB
    """

    DB_AUTH_STRING = "dbname='covid' user='covid_user' host='host_or_ip.com' port='5432' password='covid_pass'"

    def __init__(self, file_name):
        if not os.path.isfile(file_name):
            raise FileNotFoundError(f"Cannot find the specified file: '{file_name}'")
        self.df_vcf = pd.read_csv(file_name, sep=',', skiprows=0)
        self.df_vcf_sorted = None
        self.file_name = file_name
        self.full_path = (os.path.dirname(self.file_name) + '/') if os.path.dirname(self.file_name) != '' else ''

    def sort_relevant_items(self):
        print('Selecting relevant rows...')
        start = time.time()
        self.df_vcf_sorted = self.df_vcf[~self.df_vcf['pdf_json_files'].isnull() \
                & ~self.df_vcf['source_x'].str.contains('PMC') & ~self.df_vcf['url'].isnull() \
                & (self.df_vcf['abstract'].str.contains('covid', na=False, case=False) | self.df_vcf['abstract'].str.contains('sars-cov-2', na=False, case=False) \
                | self.df_vcf['title'].str.contains('covid', na=False, case=False) | self.df_vcf['title'].str.contains('sars-cov-2', na=False, case=False)) \
                ]
        print(f'Found {self.df_vcf_sorted.shape[0]} relevant rows, it took {int(time.time() - start)} sec')

    def create_stripped_archive(self, out_path):
        if self.df_vcf_sorted == None:
            self.sort_relevant_items()

        print('Copying relevant articles to the new path...')
        os.makedirs(out_path + '/' + 'document_parses', exist_ok=True)
        os.makedirs(out_path + '/' + 'document_parses/pdf_json', exist_ok=True)
        for article in self.df_vcf_sorted['pdf_json_files']:
            copyfile(self.full_path + article.split(';')[0], out_path + '/' + article.split(';')[0])
    
        print('Copying stripped metadata to the new path...')
        self.df_vcf_sorted.to_csv(out_path + '/' + 'metadata.csv')
        print('Creating stripped archive done!')

    def put_relevant_items_to_db(self):
        if self.df_vcf_sorted == None:
            self.sort_relevant_items()

        start = time.time()
        df_tmp = self.df_vcf_sorted.applymap(str)
        print('Processing META column...')
        df_tmp.replace('"', '', regex=True, inplace=True) 
        df_tmp['meta'] = df_tmp.apply(lambda row: ",".join(('"' if ',' in i else '') + i + ('"' if ',' in i else '') for i in row[:]), axis=1)
        def get_body(row):
            return open(self.full_path + row['pdf_json_files'].split(';')[0], 'r').read()

        print('Processing BODY column...')
        df_tmp['body'] = df_tmp.apply(get_body, axis=1)
        df_tmp['id'] = df_tmp.apply(lambda _: str(uuid4()), axis=1)
        df_tmp['created'] = datetime.utcnow()
        df_tmp['status'] = 'FETCHED'
        df_tmp['source'] = 'cord19_pdf'
        df_tmp['message'] = ''
        df_tmp.drop(df_tmp.columns.difference(['body', 'id', 'created', 'status', 'source', 'message', 'cord_uid', 'meta']), 1, inplace=True)
        print(f'Preparing columns took {int(time.time() - start)} sec')

        start = time.time()
        print('Connecting to DB...')
        conn = psycopg2.connect(self.DB_AUTH_STRING)
        cur = conn.cursor()

        rows = zip(df_tmp.created, df_tmp.created, df_tmp.id, df_tmp.cord_uid, df_tmp.meta, df_tmp.body, df_tmp.status, df_tmp.source, df_tmp.message)

        print('Executing query...')
        try:
            extras.execute_batch(cur, """INSERT INTO article (created, updated, id, external_id, meta, body, status, source, message) VALUES(
                                %s, %s, %s, %s, %s, %s, %s, %s, %s)  ON CONFLICT DO NOTHING""", rows, 200)
            conn.commit()
        except psycopg2.Error as e:
            print(f"DB error: {e.pgcode}")
            conn.rollback()

        cur.close()
        conn.close()
        print(f'DB operations took {int(time.time() - start)} sec')
        print('Done!')

    def get_metadata_info(self):
        print('Total rows: ', self.df_vcf.shape[0])


if __name__ == '__main__':
    cord_parser_obj = CordParser('metadata.csv')
    cord_parser_obj.put_relevant_items_to_db()

