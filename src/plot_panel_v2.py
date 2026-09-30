import os
import sys
import time
import traceback
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

print("=" * 60)
print("figures_light.py starting")
print("python:", sys.executable)
print("script:", os.path.abspath(__file__))
print("=" * 60)

plt.style.use('default')
mpl.rcParams.update({
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "savefig.facecolor": "white",
    "savefig.dpi":       300,
    "text.color":        "black",
    "axes.labelcolor":   "black",
    "axes.edgecolor":    "black",
    "xtick.color":       "black",
    "ytick.color":       "black",
    "axes.titlecolor":   "black",
    "font.size":         9,
    "legend.framealpha": 1.0,
    "legend.facecolor":  "white",
    "legend.edgecolor":  "0.7",
})

PANEL = r"D:\GENOMEINDIA\outputs\multi_gene\full_panel.csv"
STATS = r"D:\GENOMEINDIA\outputs\stats\statistical_comparison_v2.csv"
FIG   = r"D:\GENOMEINDIA\figures"
SUFFIX = "_light"

SPINE = "#999999"
NONSIG = "#A6A6A6"

POP_LABELS = {
    "india_af": "India\n(GenomeIndia)", "gnomad_sas_af": "South Asian",
    "gnomad_eas_af": "East Asian", "gnomad_afr_af": "African",
    "gnomad_nfe_af": "European", "gnomad_amr_af": "Adm. American",
    "gnomad_mid_af": "Middle Eastern",
}
POP_COLORS = {
    "india_af": "#D35400", "gnomad_sas_af": "#D68910", "gnomad_eas_af": "#2471A3",
    "gnomad_afr_af": "#1E8449", "gnomad_nfe_af": "#76448A",
    "gnomad_amr_af": "#117A65", "gnomad_mid_af": "#A93226",
}
CAT_COLORS = {"disease": "#2471A3", "pgx": "#C0392B", "control": "#4D5656"}
mkdir src outputs figuresCAT_LABELS = {"disease": "Disease susceptibility", "pgx": "Pharmacogenomic",
              "control": "Positive control"}


def outpath(stem):
    return os.path.join(FIG, f"{stem}{SUFFIX}.png")


def save(fig_path, **kw):
    plt.savefig(fig_path, bbox_inches="tight", **kw)
    plt.close()
    size = os.path.getsize(fig_path)
    when = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(fig_path)))
    print(f"  Saved {os.path.abspath(fig_path)}  ({size:,} bytes, written {when})")


def legend_patches():
    return [mpatches.Patch(color=c, label=CAT_LABELS[k]) for k, c in CAT_COLORS.items()]


def style_axes(ax):
    for s in ax.spines.values():
        s.set_color(SPINE)


def load():
    df = pd.read_csv(PANEL)
    print(f"  full_panel.csv columns: {list(df.columns)}")
    dup = df["gene"].duplicated(keep=False)
    df["label"] = np.where(dup, df["gene"] + "\n" + df["rsid"], df["gene"])
    df["short"] = np.where(
        dup,
        df["gene"] + " " + df["phenotype"].str.extract(r"\((\*\d+\w*)\)")[0].fillna(""),
        df["gene"])
    df = df.sort_values(["category", "deviation"]).reset_index(drop=True)
    return df


def pops_present(df):
    return [c for c in POP_LABELS if c in df.columns and df[c].notna().any()]


def bar_grid(df):
    pops = pops_present(df)
    n = len(df)
    ncols, nrows = 4, int(np.ceil(n / 4))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for i, row in df.iterrows():
        ax = axes[i]
        vals = [row[c] if pd.notna(row[c]) else 0.0 for c in pops]
        bars = ax.bar([POP_LABELS[c] for c in pops], vals,
                      color=[POP_COLORS[c] for c in pops],
                      edgecolor="black", linewidth=0.4)
        bars[pops.index("india_af")].set_linewidth(2.5)

        top = max(vals) if max(vals) > 0 else 0.05
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + top * 0.02,
                    f"{v:.3f}", ha="center", va="bottom", color="black", fontsize=6.5)

        ax.set_title(f"{row['label']}\n{row['phenotype']}",
                     color=CAT_COLORS[row["category"]], fontsize=9.5,
                     fontweight="bold", pad=9)
        ax.set_ylabel("Risk allele frequency", fontsize=8.5)
        ax.tick_params(axis="x", rotation=55, labelsize=6.5)
        ax.tick_params(axis="y", labelsize=7.5)
        style_axes(ax)
        ax.set_ylim(0, top * 1.22)

    for j in range(n, len(axes)):
        axes[j].axis("off")

    fig.legend(handles=legend_patches(), loc="lower right", fontsize=10,
               bbox_to_anchor=(0.98, 0.02))
    plt.tight_layout()
    save(outpath("panel_bars"), dpi=200)


def heatmap(df):
    pops = pops_present(df)
    m = df[pops].astype(float).values
    ylabels = [f"{r['short']}\n({r['phenotype'][:28]})" for _, r in df.iterrows()]

    fig, ax = plt.subplots(figsize=(10, 0.72 * len(df) + 3))
    im = ax.imshow(m, cmap="YlOrRd", aspect="auto")

    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels([POP_LABELS[c].replace("\n", " ") for c in pops],
                       rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=7.5)
    for tick, cat in zip(ax.get_yticklabels(), df["category"]):
        tick.set_color(CAT_COLORS[cat])

    vmax = np.nanmax(m)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=7.5,
                    color="white" if v > vmax * 0.65 else "black")

    ax.add_patch(plt.Rectangle((pops.index("india_af") - 0.5, -0.5), 1, len(ylabels),
                               fill=False, edgecolor="black", linewidth=2.5))
    cb = plt.colorbar(im, ax=ax, label="Risk allele frequency")
    cb.ax.yaxis.label.set_color("black")
    cb.ax.tick_params(colors="black")
    cb.outline.set_edgecolor(SPINE)
    ax.set_title("Allele frequencies across global populations\n"
                 "(row labels coloured by locus category)",
                 fontsize=12.5, fontweight="bold", pad=14)
    plt.tight_layout()
    save(outpath("panel_heatmap"))


def deviation(df):
    d = df.sort_values("deviation").copy()
    fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(d) + 3))

    bars = ax.barh(d["short"], d["deviation"],
                   color=[CAT_COLORS[c] for c in d["category"]],
                   edgecolor="black", linewidth=0.4)

    span = max(abs(d["deviation"].min()), abs(d["deviation"].max()))
    ax.set_xlim(d["deviation"].min() - span * 0.22, d["deviation"].max() + span * 0.22)
    for b, v in zip(bars, d["deviation"]):
        off = span * 0.03
        ax.text(v + off if v >= 0 else v - off, b.get_y() + b.get_height() / 2,
                f"{v:+.4f}", va="center", ha="left" if v >= 0 else "right",
                color="black", fontsize=8.5)

    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("India AF minus gnomAD global pooled AF", fontsize=10.5)
    ax.set_title("Deviation from the gnomAD global pooled frequency\n"
                 "coloured by locus category", fontsize=12.5, fontweight="bold", pad=14)
    ax.tick_params(labelsize=9)
    style_axes(ax)
    ax.legend(handles=legend_patches(), fontsize=9, loc="lower right")
    plt.tight_layout()
    save(outpath("panel_deviation"))


def forest(stats_df):
    """
    Odds ratios with 95% CI on a log scale, one row per locus-population
    comparison, grouped by category. Non-significant rows are drawn in grey
    rather than at low alpha, which nearly vanishes on a white background.
    """
    s = stats_df.copy()
    dup = s.groupby("gene")["rsid"].transform("nunique") > 1
    s["short"] = np.where(dup, s["gene"] + " " + s["rsid"], s["gene"])
    s["ylab"] = s["short"] + " — " + s["population"]

    cat_order = {"control": 0, "disease": 1, "pgx": 2}
    s["_o"] = s["category"].map(cat_order)
    s = s.sort_values(["_o", "short", "odds_ratio"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 0.23 * len(s) + 3))

    for i, r in s.iterrows():
        sig = bool(r["significant_fdr"])
        col = CAT_COLORS[r["category"]] if sig else NONSIG
        ax.plot([r["or_lower"], r["or_upper"]], [i, i],
                color=col, lw=1.4 if sig else 1.0, zorder=2 if sig else 1)
        ax.plot(r["odds_ratio"], i, "o", color=col,
                ms=4.0 if sig else 3.0, zorder=3 if sig else 1)

    ax.axvline(1, color="black", lw=1, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(s)))
    ax.set_yticklabels(s["ylab"], fontsize=5.8)
    for tick, (cat, sig) in zip(ax.get_yticklabels(),
                                zip(s["category"], s["significant_fdr"])):
        tick.set_color(CAT_COLORS[cat] if sig else NONSIG)
    ax.set_xlabel("Odds ratio, India vs comparison population (log scale)", fontsize=10.5)
    ax.set_title("Odds ratios with 95% confidence intervals\n"
                 "points right of the dashed line are elevated in India; "
                 "grey rows are not significant after FDR correction",
                 fontsize=12, fontweight="bold", pad=14)
    ax.grid(axis="x", color="0.9", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    style_axes(ax)
    ax.legend(handles=legend_patches() + [mpatches.Patch(color=NONSIG, label="Not significant")],
              fontsize=9, loc="lower right")
    ax.set_ylim(-1, len(s))
    plt.tight_layout()
    save(outpath("panel_forest"))


def elevation_count(stats_df):
    """
    The paper's central figure: how many populations is each locus
    significantly elevated against? Makes the 6/6 and 4/6 immediately visible.
    """
    s = stats_df.copy()
    dup = s.groupby("gene")["rsid"].transform("nunique") > 1
    s["short"] = np.where(dup, s["gene"] + " " + s["rsid"], s["gene"])
    sig = s[s["significant_fdr"]]
    g = (sig.assign(hi=sig["direction"].eq("higher"))
            .groupby(["category", "short"])["hi"].sum().reset_index())
    cat_order = {"control": 0, "disease": 1, "pgx": 2}
    g["_o"] = g["category"].map(cat_order)
    g = g.sort_values(["_o", "hi"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9, 0.5 * len(g) + 3))
    bars = ax.barh(g["short"], g["hi"],
                   color=[CAT_COLORS[c] for c in g["category"]],
                   edgecolor="black", linewidth=0.4)
    for b, v in zip(bars, g["hi"]):
        ax.text(v + 0.08, b.get_y() + b.get_height() / 2, f"{int(v)}/6",
                va="center", color="black", fontsize=9)

    ax.axvline(4, color="black", lw=1, ls=":")
    ax.text(4.05, -0.7, "majority rule (4 of 6)", fontsize=8, color="black")

    ax.set_xlim(0, 7)
    ax.set_xticks(range(7))
    ax.set_xlabel("Populations where India is significantly higher (of 6)", fontsize=10.5)
    ax.set_title("Elevation breadth by locus\n"
                 "no disease locus exceeds 3 of 6; two pharmacogenomic loci do",
                 fontsize=12.5, fontweight="bold", pad=14)
    ax.tick_params(labelsize=9)
    for tick, cat in zip(ax.get_yticklabels(), g["category"]):
        tick.set_color(CAT_COLORS[cat])
    ax.grid(axis="x", color="0.9", lw=0.6)
    ax.set_axisbelow(True)
    style_axes(ax)
    ax.legend(handles=legend_patches(), fontsize=9, loc="lower right")
    plt.tight_layout()
    save(outpath("panel_elevation"))


def main():
    for p in (PANEL, STATS):
        if not os.path.exists(p):
            sys.exit(f"MISSING INPUT: {p}")
    os.makedirs(FIG, exist_ok=True)
    print(f"Output directory: {os.path.abspath(FIG)}\n")

    print("Loading panel ...")
    df = load()
    print(f"  Loci: {len(df)}")
    print("Loading stats ...")
    stats_df = pd.read_csv(STATS)
    print(f"  Comparisons: {len(stats_df)}\n")

    jobs = [
        ("panel_bars",      lambda: bar_grid(df)),
        ("panel_heatmap",   lambda: heatmap(df)),
        ("panel_deviation", lambda: deviation(df)),
        ("panel_forest",    lambda: forest(stats_df)),
        ("panel_elevation", lambda: elevation_count(stats_df)),
    ]
    ok, failed = 0, []
    for name, fn in jobs:
        print(f"Building {name}{SUFFIX}.png ...")
        try:
            fn()
            ok += 1
        except Exception:
            failed.append(name)
            print(f"  FAILED: {name}")
            traceback.print_exc()
            plt.close("all")

    print("\n" + "=" * 60)
    print(f"Done. {ok} of {len(jobs)} figures written to {os.path.abspath(FIG)}")
    if failed:
        print("Failed:", ", ".join(failed))
    print("=" * 60)


if __name__ == "__main__":
    main()