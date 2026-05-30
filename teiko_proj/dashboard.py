import io
import base64
 
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
 
import dash
from dash import dcc, html, dash_table, Input, Output, State
import dash_table
 
from analysis import (
    get_frequency_table,
    get_responder_data,
    run_statistics,
    get_baseline_subset,
)
 
RESP_COLORS   = {"yes": "#27ae60", "no": "#e74c3c"}
POP_COLORS    = {
    "b_cell":       "#3498db",
    "cd8_t_cell":   "#9b59b6",
    "cd4_t_cell":   "#e67e22",
    "nk_cell":      "#1abc9c",
    "monocyte":     "#e74c3c",
}
BG            = "#f0f2f5"
CARD_BG       = "#ffffff"
HEADER_COLOR  = "#1a2c45"
ACCENT        = "#2980b9"
 
freq_df   = get_frequency_table()
resp_df   = get_responder_data()
stats_df  = run_statistics()
subset    = get_baseline_subset()
 
POPULATIONS = sorted(freq_df["population"].unique())
SAMPLES     = sorted(freq_df["sample"].unique())
 
 
def card(children, style=None):
    base = {
        "background": CARD_BG,
        "borderRadius": "10px",
        "padding": "24px",
        "marginBottom": "20px",
        "boxShadow": "0 2px 8px rgba(0,0,0,0.08)",
    }
    if style:
        base.update(style)
    return html.Div(children, style=base)
 
 
def section_title(text):
    return html.H3(text, style={
        "color": HEADER_COLOR,
        "borderBottom": f"3px solid {ACCENT}",
        "paddingBottom": "8px",
        "marginTop": "0",
    })
 
 
 
def fig_frequency_bar(sample_ids=None):
    df = freq_df.copy()
    if sample_ids:
        df = df[df["sample"].isin(sample_ids)]
    fig = px.bar(
        df, x="sample", y="percentage", color="population",
        color_discrete_map=POP_COLORS,
        labels={"percentage": "Relative Frequency (%)", "sample": "Sample", "population": "Population"},
        title="Cell Population Frequencies per Sample",
        barmode="stack",
    )
    fig.update_layout(
        plot_bgcolor=CARD_BG, paper_bgcolor=CARD_BG,
        legend_title="Population",
        xaxis_tickangle=-45,
        margin=dict(b=120),
    )
    return fig
 
 
def fig_boxplots():
    pops = POPULATIONS
    fig  = make_subplots(rows=1, cols=len(pops),
                         subplot_titles=[p.replace("_"," ").title() for p in pops])
 
    stat_map = stats_df.set_index("population")
 
    for i, pop in enumerate(pops, start=1):
        for resp, label, color in [("yes","Responder","#27ae60"), ("no","Non-Resp","#e74c3c")]:
            vals = resp_df[(resp_df["population"]==pop) & (resp_df["response"]==resp)]["percentage"].values
            fig.add_trace(
                go.Box(
                    y=vals, name=label, marker_color=color,
                    boxmean=True, showlegend=(i==1),
                    legendgroup=label,
                    hovertemplate=f"<b>{label}</b><br>%{{y:.2f}}%<extra></extra>",
                ),
                row=1, col=i,
            )
            
            fig.add_trace(
                go.Scatter(
                    x=np.random.uniform(-0.2, 0.2, len(vals)),
                    y=vals,
                    mode="markers",
                    marker=dict(color=color, size=4, opacity=0.4),
                    showlegend=False,
                    hoverinfo="skip",
                ),
                row=1, col=i,
            )
 
        pval = stat_map.loc[pop, "p_value"]
        sig  = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        fig.add_annotation(
            text=f"p={pval:.3f} {sig}", xref="paper",
            x=(i - 0.5) / len(pops), y=1.07, yref="paper",
            showarrow=False, font=dict(size=11, color="navy"),
        )
 
    fig.update_layout(
        title="Responders vs Non-Responders by Cell Population<br><sup>Melanoma · miraclib · PBMC</sup>",
        plot_bgcolor=CARD_BG, paper_bgcolor=CARD_BG,
        height=480,
        legend=dict(orientation="h", y=-0.18),
    )
    return fig
 
 
def fig_stats_table():
    df = stats_df.copy()
    df["median_responders"]     = df["median_responders"].round(3)
    df["median_non_responders"] = df["median_non_responders"].round(3)
    df["p_value"]               = df["p_value"].round(5)
    df["mann_whitney_u"]        = df["mann_whitney_u"].astype(int)
    df["significant (p<0.05)"]  = df["significant (p<0.05)"].map({True: "✅ Yes", False: "No"})
    df["population"]            = df["population"].str.replace("_", " ").str.title()
    return df
 
 
def fig_subset_sunburst():
    sub_df = subset["samples"].copy()
    sub_df["sex_label"] = sub_df["sex"].map({"M": "Male", "F": "Female"})
    sub_df["resp_label"] = sub_df["response"].map({"yes": "Responder", "no": "Non-Resp"})
    grouped = sub_df.groupby(["project","resp_label","sex_label"]).size().reset_index(name="n")
    fig = px.sunburst(
        grouped, path=["project","resp_label","sex_label"], values="n",
        title="Baseline Samples: Project → Response → Sex",
        color="project", color_discrete_sequence=px.colors.qualitative.Pastel,
    )
    fig.update_layout(paper_bgcolor=CARD_BG, height=420)
    return fig
 
 
def fig_subset_bars():
    proj = subset["samples_per_project"]
    fig  = px.bar(proj, x="project", y="n_samples",
                  color="project", text="n_samples",
                  title="Baseline Samples per Project",
                  color_discrete_sequence=px.colors.qualitative.Set2)
    fig.update_traces(textposition="outside")
    fig.update_layout(showlegend=False, plot_bgcolor=CARD_BG, paper_bgcolor=CARD_BG)
    return fig
 
 
 
app = dash.Dash(__name__, title="Loblaw Bio · Clinical Trial Dashboard")
app.config.suppress_callback_exceptions = True
 
app.layout = html.Div(style={"backgroundColor": BG, "minHeight": "100vh",
                               "fontFamily": "Inter, Segoe UI, sans-serif"}, children=[
 
    
    html.Div(style={
        "background": f"linear-gradient(135deg, {HEADER_COLOR} 0%, #2c4a70 100%)",
        "padding": "28px 40px", "color": "white", "marginBottom": "28px",
    }, children=[
        html.H1("🔬 Loblaw Bio · Clinical Trial Dashboard",
                style={"margin": "0", "fontSize": "26px", "fontWeight": "700"}),
        html.P("Immune cell population analysis · miraclib drug candidate",
               style={"margin": "6px 0 0", "opacity": "0.75", "fontSize": "14px"}),
    ]),
 
    html.Div(style={"maxWidth": "1400px", "margin": "0 auto", "padding": "0 24px"}, children=[
 
        
        dcc.Tabs(id="tabs", value="tab-overview", style={"marginBottom": "20px"},
                 colors={"border": BG, "primary": ACCENT, "background": CARD_BG},
                 children=[
                     dcc.Tab(label="📊 Data Overview",       value="tab-overview"),
                     dcc.Tab(label="📈 Statistical Analysis", value="tab-stats"),
                     dcc.Tab(label="🔍 Subset Analysis",      value="tab-subset"),
                 ]),
 
        html.Div(id="tab-content"),
    ]),
])
 
 
 
@app.callback(Output("tab-content", "children"), Input("tabs", "value"))
def render_tab(tab):
 
    
    if tab == "tab-overview":
        return html.Div([
            card([
                section_title("Part 2 · Cell Population Frequency Overview"),
                html.P("Select samples to filter the chart and frequency table below.",
                       style={"color": "#666"}),
                html.Div([
                    html.Label("Filter by sample:", style={"fontWeight": "600", "marginRight": "10px"}),
                    dcc.Dropdown(
                        id="sample-filter",
                        options=[{"label": s, "value": s} for s in SAMPLES],
                        multi=True, placeholder="All samples",
                        style={"flexGrow": "1"},
                    ),
                ], style={"display": "flex", "alignItems": "center", "marginBottom": "16px"}),
                dcc.Graph(id="freq-bar"),
            ]),
            card([
                section_title("Frequency Table"),
                html.Div([
                    html.Label("Filter by population:", style={"fontWeight": "600", "marginRight": "10px"}),
                    dcc.Dropdown(
                        id="pop-filter",
                        options=[{"label": p.replace("_"," ").title(), "value": p} for p in POPULATIONS],
                        multi=True, placeholder="All populations",
                        style={"width": "350px"},
                    ),
                ], style={"display": "flex", "alignItems": "center", "marginBottom": "16px"}),
                html.Div(id="freq-table"),
            ]),
        ])
 
    
    elif tab == "tab-stats":
        sig_pops = stats_df[stats_df["significant (p<0.05)"]]["population"].tolist()
        sig_text = (
            f"Statistically significant differences (p < 0.05) were found in: "
            f"{', '.join(p.replace('_',' ').title() for p in sig_pops)}."
            if sig_pops else
            "No populations reached statistical significance at p < 0.05."
        )
        return html.Div([
            card([
                section_title("Part 3 · Responders vs Non-Responders"),
                html.Div(style={
                    "background": "#eaf4fb", "borderLeft": f"4px solid {ACCENT}",
                    "padding": "14px 18px", "borderRadius": "4px", "marginBottom": "20px",
                }, children=[
                    html.P([
                        html.Strong("Cohort: "),
                        "Melanoma patients · miraclib treatment · PBMC samples only"
                    ], style={"margin": "0 0 6px"}),
                    html.P([html.Strong("Statistical test: "), "Mann-Whitney U (two-sided, non-parametric)"],
                           style={"margin": "0 0 6px"}),
                    html.P([html.Strong("Finding: "), sig_text], style={"margin": "0"}),
                ]),
                dcc.Graph(figure=fig_boxplots()),
            ]),
            card([
                section_title("Statistical Test Results"),
                dash_table.DataTable(
                    data=fig_stats_table().to_dict("records"),
                    columns=[{"name": c.replace("_", " ").title(), "id": c}
                             for c in fig_stats_table().columns],
                    style_table={"overflowX": "auto"},
                    style_header={
                        "backgroundColor": HEADER_COLOR, "color": "white",
                        "fontWeight": "600", "textAlign": "center",
                    },
                    style_cell={"textAlign": "center", "padding": "10px", "fontSize": "13px"},
                    style_data_conditional=[
                        {
                            "if": {"filter_query": '{significant (p<0.05)} = "✅ Yes"'},
                            "backgroundColor": "#d5f5e3", "fontWeight": "600",
                        },
                    ],
                ),
            ]),
        ])
 
    
    elif tab == "tab-subset":
        sub_df = subset["samples"]
        
        summary_cards = html.Div(style={"display": "grid",
                                        "gridTemplateColumns": "repeat(4, 1fr)",
                                        "gap": "16px", "marginBottom": "20px"}, children=[
            _metric_card("Total Samples",  str(subset["total_samples"]),  "🧪"),
            _metric_card("Total Subjects", str(subset["total_subjects"]), "👤"),
            _metric_card("Responders",
                         str(subset["response_counts"].set_index("response_status").loc["Responders","n_subjects"]),
                         "✅"),
            _metric_card("Non-Responders",
                         str(subset["response_counts"].set_index("response_status").loc["Non-Responders","n_subjects"]),
                         "❌"),
        ])
 
        return html.Div([
            card([
                section_title("Part 4 · Baseline Subset Analysis"),
                html.Div(style={
                    "background": "#fef9e7", "borderLeft": "4px solid #f39c12",
                    "padding": "14px 18px", "borderRadius": "4px", "marginBottom": "20px",
                }, children=[
                    html.P([html.Strong("Filter: "), "Melanoma · miraclib · PBMC · time_from_treatment_start = 0"],
                           style={"margin": "0"}),
                ]),
                summary_cards,
            ]),
            html.Div(style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "20px"}, children=[
                card([dcc.Graph(figure=fig_subset_bars())]),
                card([dcc.Graph(figure=fig_subset_sunburst())]),
            ]),
            card([
                section_title("Response & Sex Breakdown"),
                html.Div(style={"display": "grid", "gridTemplateColumns": "1fr 1fr", "gap": "24px"}, children=[
                    html.Div([
                        html.H4("Response Status", style={"marginTop": "0"}),
                        dash_table.DataTable(
                            data=subset["response_counts"].to_dict("records"),
                            columns=[{"name": c.replace("_"," ").title(), "id": c}
                                     for c in subset["response_counts"].columns],
                            style_header={"backgroundColor": HEADER_COLOR, "color": "white", "fontWeight": "600"},
                            style_cell={"textAlign": "center", "padding": "10px"},
                        ),
                    ]),
                    html.Div([
                        html.H4("Sex Distribution", style={"marginTop": "0"}),
                        dash_table.DataTable(
                            data=subset["sex_counts"].to_dict("records"),
                            columns=[{"name": c.replace("_"," ").title(), "id": c}
                                     for c in subset["sex_counts"].columns],
                            style_header={"backgroundColor": HEADER_COLOR, "color": "white", "fontWeight": "600"},
                            style_cell={"textAlign": "center", "padding": "10px"},
                        ),
                    ]),
                ]),
            ]),
            card([
                section_title("Raw Baseline Sample Records"),
                dash_table.DataTable(
                    data=sub_df.to_dict("records"),
                    columns=[{"name": c.replace("_"," ").title(), "id": c} for c in sub_df.columns],
                    page_size=15,
                    sort_action="native",
                    filter_action="native",
                    style_table={"overflowX": "auto"},
                    style_header={"backgroundColor": HEADER_COLOR, "color": "white", "fontWeight": "600"},
                    style_cell={"textAlign": "center", "padding": "8px", "fontSize": "12px"},
                    style_data_conditional=[
                        {"if": {"row_index": "odd"}, "backgroundColor": "#f7f9fc"},
                    ],
                ),
            ]),
        ])
 
 
def _metric_card(label, value, icon):
    return html.Div(style={
        "background": CARD_BG, "borderRadius": "10px",
        "padding": "20px", "textAlign": "center",
        "boxShadow": "0 2px 8px rgba(0,0,0,0.07)",
    }, children=[
        html.Div(icon, style={"fontSize": "28px", "marginBottom": "6px"}),
        html.Div(value, style={"fontSize": "32px", "fontWeight": "700", "color": HEADER_COLOR}),
        html.Div(label, style={"fontSize": "13px", "color": "#888", "marginTop": "4px"}),
    ])
 

 
@app.callback(
    Output("freq-bar", "figure"),
    Input("sample-filter", "value"),
)
def update_bar(sample_ids):
    return fig_frequency_bar(sample_ids if sample_ids else None)
 
 
@app.callback(
    Output("freq-table", "children"),
    [Input("sample-filter", "value"), Input("pop-filter", "value")],
)
def update_table(sample_ids, pop_ids):
    df = freq_df.copy()
    if sample_ids:
        df = df[df["sample"].isin(sample_ids)]
    if pop_ids:
        df = df[df["population"].isin(pop_ids)]
 
    df_display = df.copy()
    df_display["percentage"] = df_display["percentage"].round(2).astype(str) + " %"
    df_display["population"] = df_display["population"].str.replace("_"," ").str.title()
 
    return dash_table.DataTable(
        data=df_display.head(500).to_dict("records"),
        columns=[{"name": c.replace("_"," ").title(), "id": c} for c in df_display.columns],
        page_size=20,
        sort_action="native",
        filter_action="native",
        style_table={"overflowX": "auto"},
        style_header={"backgroundColor": HEADER_COLOR, "color": "white", "fontWeight": "600"},
        style_cell={"textAlign": "center", "padding": "8px", "fontSize": "13px"},
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#f7f9fc"},
        ],
    )
 
 
if __name__ == "__main__":
    app.run(debug=True, port=8050)