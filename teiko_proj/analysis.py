import sqlite3
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

DB_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'clinical_trial.db')

def get_connection():
    return sqlite3.connect(DB_PATH)

def get_frequency_table():
    conn = get_connection()
    query = """"
        SELECT 
            cc.sample,
            cc.population,
            cc.count,
            SUM(cc.count) OVER (PARTITION BY cc.sample) AS total_count,
            ROUND(100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample), 4) AS percentage
            FROM cell_counts cc
            ORDER BY cc.sample, cc.population
        """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df[["sample", "total_count", "population", "count", "percentage"]]

def get_responder_data():
    conn = get_connection()
    query = """
        SELECT 
            cc.sample,
            su.subject,
            su.response,
            cc.population,
            cc.count,
            SUM(cc.count) OVER (PARTITION BY cc.sample) AS total_count,
            ROUND(100.0 * cc.count / SUM(cc.count) OVER (PARTITION BY cc.sample), 4) AS percentage
        FROM cell_counts cc
        JOIN samples sa ON cc.sample = sa.sample
        JOIN subjects su ON sa.subject = su.subject
        WHERE su.condition = 'melanoma'
            AND sa.treatment = 'miraclib'
            AND su.sample_type = 'PBMC'
            AND su.response IS NOT NULL
        ORDER BY cc.sample, cc.population
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df

if __name__ == "__main__":
    freq_table = get_frequency_table()
    print(freq_table.head(10).to_string(index=False))