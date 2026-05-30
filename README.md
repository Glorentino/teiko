# Clinical Trial Dashboard

An interactive data analysis dashboard for exploring immune cell population data from a clinical trial of the drug candidate **miraclib** in melanoma patients. The pipeline ingests raw cell-count data, runs statistical comparisons between responders and non-responders, and visualizes results in a Plotly Dash web application.

---

## How to Run

All commands should be run from the **root of the repository**.

### 1. Install dependencies

```bash
make setup
```

This runs `pip install -r teiko_proj/requirements.txt`.

### 2. Run the full data pipeline

```bash
make pipeline
```

This runs two steps sequentially:
1. `load_data.py` — initializes the SQLite database and loads `teiko_proj/cell-count.csv`
2. `teiko_proj/analysis.py` — runs all analyses and saves output files to `outputs/`

Generated output files:

| File | Description |
|---|---|
| `outputs/frequency_table.csv` | Relative frequency (%) of each cell population per sample |
| `outputs/stats_results.csv` | Mann-Whitney U test results per population |
| `outputs/boxplot.png` | Box plot comparing responders vs non-responders |
| `outputs/baseline_subset.csv` | Raw baseline (t=0) sample records |
| `outputs/baseline_samples_per_project.csv` | Sample counts per project |
| `outputs/baseline_response_counts.csv` | Responder / non-responder subject counts |
| `outputs/baseline_sex_counts.csv` | Sex distribution of baseline subjects |

### 3. Start the dashboard

```bash
make dashboard
```

Opens the Dash app at **http://localhost:8050**.

> **GitHub Codespaces:** After running `make dashboard`, Codespaces will automatically forward port 8050. Click the pop-up notification or go to the **Ports** tab and open the forwarded URL.

---

## Dashboard Link

Run `make dashboard` and open: http://localhost:8050

---

## Database Schema

### Tables

#### `subjects`
| Column | Type | Description |
|---|---|---|
| subject | TEXT (PK) | Unique subject identifier |
| age | INTEGER | Subject age |
| sex | TEXT | `M` or `F` |
| condition | TEXT | Disease condition (e.g. `melanoma`, `carcinoma`) |
| response | TEXT | Drug response: `yes`, `no`, or NULL if unknown |

#### `projects`
| Column | Type | Description |
|---|---|---|
| project | TEXT (PK) | Unique project identifier |

#### `samples`
| Column | Type | Description |
|---|---|---|
| sample | TEXT (PK) | Unique sample identifier |
| subject | TEXT (FK) | References `subjects.subject` |
| project | TEXT (FK) | References `projects.project` |
| sample_type | TEXT | Biospecimen type (e.g. `PBMC`) |
| treatment | TEXT | Treatment administered (e.g. `miraclib`) |
| time_from_treatment_start | INTEGER | Days since treatment start; 0 = pre-treatment baseline |

#### `cell_counts`
| Column | Type | Description |
|---|---|---|
| id | INTEGER (PK) | Auto-incremented row ID |
| sample | TEXT (FK) | References `samples.sample` |
| population | TEXT | Cell population name (`b_cell`, `cd8_t_cell`, `cd4_t_cell`, `nk_cell`, `monocyte`) |
| count | INTEGER | Raw cell count for this population in this sample |

### Design Rationale

The flat CSV is normalized into four tables to eliminate redundancy and make querying flexible:

- **subjects** and **projects** are reference tables. Storing them separately means subject demographics and project metadata are never duplicated across rows.
- **samples** is the join point between a subject, a project, and a collection event. The `time_from_treatment_start` column on `samples` (rather than on `subjects`) correctly models that the same subject can have multiple samples at different time points.
- **cell_counts** is the fact table. Storing one row per (sample, population) pair — rather than one wide row with five population columns — makes it straightforward to add new populations without a schema migration, and allows efficient filtering by population via a simple `WHERE` clause.
- Relative frequencies are **computed at query time** using a SQL window function (`SUM(...) OVER (PARTITION BY sample)`) rather than stored. This avoids stale derived values if raw counts are ever corrected.

### Scalability

| Concern | Approach |
|---|---|
| Hundreds of projects / thousands of samples | The normalized schema already handles this; add indexes on `samples(subject)`, `samples(project)`, and `cell_counts(sample, population)` to keep joins fast |
| New cell populations | Insert new rows into `cell_counts` — no schema change needed |
| New analytics (longitudinal, survival, clustering) | All time-series queries use `time_from_treatment_start` already present on `samples`; survival / clustering layers can join against `subjects` without touching the fact table |
| Concurrent reads | Swap SQLite for PostgreSQL; the query layer (`analysis.py`) only uses standard SQL and `pandas.read_sql_query`, so the connection string is the only change |
| Large datasets | Partition `cell_counts` by `sample` (or by project via a join) and use columnar storage (e.g. Parquet + DuckDB) for analytical queries; the Python layer is unchanged |

---

## Code Structure

```
teiko/
├── load_data.py              # Part 1 — ETL: CSV → SQLite
├── clinical_trial.db         # SQLite database (generated by pipeline)
├── Makefile                  # setup / pipeline / dashboard targets
├── teiko_proj/
│   ├── cell-count.csv        # Raw input data
│   ├── requirements.txt      # Python dependencies
│   ├── analysis.py           # Parts 2–4 — queries, statistics, output files
│   └── dashboard.py          # Interactive Dash web application
└── outputs/                  # Generated by `make pipeline`
    ├── frequency_table.csv
    ├── stats_results.csv
    ├── boxplot.png
    ├── baseline_subset.csv
    ├── baseline_samples_per_project.csv
    ├── baseline_response_counts.csv
    └── baseline_sex_counts.csv
```

### `load_data.py` (Part 1)
Entry point for the ETL step. Reads `teiko_proj/cell-count.csv`, normalizes it into the four-table schema, and writes `clinical_trial.db` to the project root. Reference tables use `INSERT OR IGNORE` so re-runs are idempotent; `cell_counts` is fully wiped and reloaded to avoid duplicates.

### `teiko_proj/analysis.py` (Parts 2–4)
Stateless query and computation layer. Each function opens its own database connection and returns a pandas DataFrame or dict, making them independently testable and safe to call from the dashboard. When run as `__main__` (via `make pipeline`), it saves all outputs to `outputs/`.

- `get_frequency_table()` — relative frequency per sample/population via SQL window function
- `get_responder_data()` — cohort-filtered data (melanoma · miraclib · PBMC) joined with response status
- `run_statistics()` — two-sided Mann-Whitney U test per population, sorted by p-value
- `make_boxplot()` — matplotlib box plot with jitter and significance annotations
- `get_baseline_subset()` — pre-treatment (t=0) cohort with response and sex breakdowns

### `teiko_proj/dashboard.py` (Parts 2–4 — interactive)
Plotly Dash application with three tabs:

| Tab | Content |
|---|---|
| Data Overview | Stacked bar chart and frequency table, both filterable by sample and/or population |
| Statistical Analysis | Interactive box plots (Plotly), cohort summary callout, and statistics DataTable |
| Subset Analysis | Metric cards (total samples, subjects, responders, non-responders), bar chart, sunburst chart, and raw record table |

The app is structured so all heavy computation happens at startup (module level), and callbacks only filter already-loaded DataFrames — keeping interactions fast regardless of dataset size.
