import sqlite3
import pandas as pd
from scipy import stats
import numpy as np
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

def run_statistics():
    df = get_responder_data()
    populations = df['population'].unique()
    results = []
    for pop in populations:
        sub = df[df['population'] == pop]
        responders = sub[sub["response"] == "yes"]["percentage"].values
        non_responders = sub[sub["response"] == "no"]["percentage"].values

        stat,pval = stats.mannwhitneyu(responders, non_responders, alternative='two-sided')

        results.append({
            "population": pop,
            "n_responders": len(responders),
            "n_non_responders": len(non_responders),
            "median_responders": np.median(responders),
            "median_non_responders": np.median(non_responders),
            "mannwhitneyu_stat": stat,
            "p_value": pval,
            "significant (p < 0.05)": pval < 0.05,
        })
    results_df = pd.DataFrame(results).sort_values("p_value")
    return results_df

def make_boxplot():
    df = get_responder_data()
    populations = sorted(df['population'].unique())
    stats_df = run_statistics()

    COLORS = {"yes": "#2ecc71", "no": "#e74c3c"}
    LABELS = {"yes": "Responders", "no": "Non-Responders"}

    fig, axes = plt.subplots(1, len(populations), figsize=(16, 6), sharey=False)
    fig.patch.set_facecolor("#f8f9fa")

    for ax, pop in zip(axes, populations):
        data_yes = df[(df["population"] == pop) & (df["response"] == "yes")]["percentage"].values
        data_no = df[(df["population"] == pop) & (df["response"] == "no")]["percentage"].values

        bp = ax.boxplot(
            [data_yes, data_no],
            patch_artist=True,
            widths=0.6,
            medianprops=dict(color="black", linewidth= 2),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            flierprops=dict(marker='o', color='black', markersize=4, alpha=0.5)
        )
        bp["boxes"][0].set_facecolor(COLORS["yes"])
        bp["boxes"][1].set_facecolor(COLORS["no"])
        bp["boxes"][0].set_alpha(0.75)
        bp["boxes"][1].set_alpha(0.75)

        for i, data in enumerate([data_yes, data_no], start=1):
            jitter = np.random.uniform(-0.15, 0.15, size=len(data))
            color = COLORS["yes"] if i == 1 else COLORS["no"]
            ax.scatter(np.full(len(data), i) + jitter, data, 
                       alpha=0.5, s=20, color=color, zorder=3)
            
        row = stats_df[stats_df["population"] == pop].iloc[0]
        pval = row["p_value"]
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        ymax = max(data_yes.max() if len(data_yes) else 0,
                   data_no.max() if len(data_no) else 0)
        ax.annotate(f"p={pval:.3e} {sig}", xy=(1.5, ymax * 1.05), 
                    ha="center", fontsize=8.5, color="navy")
        
        ax.set_title(pop.replace("_", " ").title(), fontsize=11, fontweight="bold")
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Responders", "Non-\nResponders"], fontsize=9)
        ax.ylabel("Relative Frequency (%)" if pop == populations[0] else "")
        ax.grid(axis="y", linestyle="--", alpha=0.4)
    
    legend_patches = [
        mpatches.Patch(color=COLORS["yes"], alpha=0.75, label="Responders"),
        mpatches.Patch(color=COLORS["no"], alpha=0.75, label="Non-Responders")
    ]
    fig.legend(handles=legend_patches, loc="upper right", fontsize=10)
    fig.suptitle(
        "Cell Population Frequencies: Responders vs Non-Responders\n"
        "(Melanoma patients . miraclib . PBMC)",
        fontsize=13, fontweight="bold", y=1.01
    )
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    print("==== Part 2: Frequency Table (first 10 rows) ====")
    freq_table = get_frequency_table()
    print(freq_table.head(10).to_string(index=False))

    print("==== Part 3: Statistical Analysis ====")
    st = run_statistics()
    print(st.to_string(index=False))