import sqlite3
import pandas as pd
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(_ROOT, 'clinical_trial.db')
CSV_PATH = os.path.join(_ROOT, 'teiko_proj', 'cell-count.csv')

def init_db(conn):
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS subjects (
                    subject TEXT PRIMARY KEY,
                    age INTEGER,
                    sex TEXT,
                    condition TEXT,
                    response TEXT
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS projects (
                    project TEXT PRIMARY KEY
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS samples (
                sample TEXT PRIMARY KEY,
                subject TEXT NOT NULL,
                project TEXT NOT NULL,
                sample_type TEXT,
                treatment TEXT,
                time_from_treatment_start INTEGER,
                FOREIGN KEY (subject) REFERENCES subjects(subject),
                FOREIGN KEY (project) REFERENCES projects(project)
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS cell_counts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            FOREIGN KEY (sample) REFERENCES samples(sample)
        )
    """)

    conn.commit()
    print("Database initialized successfully.")

def load_csv(conn, csv_path):
    df = pd.read_csv(csv_path)
    cursor = conn.cursor()

    populations = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

    subject_cols = ["subject", "age", "sex", "condition", "response"]
    subjects = df[subject_cols].drop_duplicates(subset=["subject"])
    for _, row in subjects.iterrows():
        cursor.execute("""
            INSERT OR IGNORE INTO subjects (subject, age, sex, condition, response)
            VALUES (?, ?, ?, ?, ?)
        """, (row["subject"], row["age"], row["sex"], row["condition"], row["response"] if pd.notna(row["response"]) else None))

    for project in df["project"].unique():
        cursor.execute("INSERT OR IGNORE INTO projects (project) VALUES (?)", (project,))

    sample_cols = ["sample", "subject", "project", "sample_type", "treatment", "time_from_treatment_start"]
    samples = df[sample_cols].drop_duplicates(subset=["sample"])

    for _, row in samples.iterrows():
        cursor.execute("""
            INSERT OR IGNORE INTO samples
                (sample, subject, project, sample_type, treatment, time_from_treatment_start)
            VALUES (?, ?, ?, ?, ?, ?)
            """, (row["sample"], row["subject"], row["project"], row["sample_type"],
                row["treatment"], row["time_from_treatment_start"]))

    cursor.execute("DELETE FROM cell_counts")
    records = []
    for _, row in df.iterrows():
        for pop in populations:
            if pd.notna(row[pop]):
                records.append((row["sample"], pop, int(row[pop])))
    cursor.executemany("""
        INSERT INTO cell_counts (sample, population, count)
        VALUES (?, ?, ?)
    """, records)

    conn.commit()
    print(f"Loaded {len(df)} CSV rows -> "
          f"{len(subjects)} subjects, {len(samples)} samples, {len(records)} cell counts.")


if __name__ == "__main__":
    conn = sqlite3.connect(DB_PATH)
    init_db(conn)
    load_csv(conn, CSV_PATH)
    conn.close()
    print(f"Data loading complete. Database saved to {DB_PATH}.")
