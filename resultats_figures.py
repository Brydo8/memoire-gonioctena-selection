#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reproduit l'ensemble des figures et tableaux de la partie Résultats du mémoire
(sélection positive entre Gonioctena quinquepunctata et G. intermedia).

Entrées : les deux tableaux de sortie du pipeline
    - kaks_results.tsv      (colonnes : gene, Ka, Ks, Ka_Ks, p_fisher, p_fisher_adj, status, ...)
    - mk_results_mkado.tsv  (colonnes : gene, Dn, Ds, Pn, Ps, alpha, DoS, p_value, p_value_adjusted, ...)

Usage :
    python resultats_figures.py kaks_results.tsv mk_results_mkado.tsv [dossier_sortie]

Sorties (PNG 300 dpi) + deux CSV (15 gènes MK ; 15 meilleurs Ka/Ks).
Dépendances : pandas, numpy, matplotlib.
"""
import sys, os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

RED, GREEN, GREY, NAVY, BLUE, GOLD = "#C00000", "#548235", "#B4B4B4", "#1F3864", "#2E75B6", "#BF8F00"

def main(kaks_path, mk_path, out="."):
    os.makedirs(out, exist_ok=True)
    k = pd.read_csv(kaks_path, sep="\t"); m = pd.read_csv(mk_path, sep="\t")
    for c in ["Ka","Ks","Ka_Ks","p_fisher","p_fisher_adj"]: k[c] = pd.to_numeric(k[c], errors="coerce")
    for c in ["Dn","Ds","Pn","Ps","alpha","DoS","p_value","p_value_adjusted"]: m[c] = pd.to_numeric(m[c], errors="coerce")

    # ----- critères de sélection (voir 3.10) -----
    kaks_pos = k[(k["Ka_Ks"] > 1) & np.isfinite(k["Ka_Ks"])]
    mk_sig   = m[(m["alpha"] > 0) & (m["p_value"] < 0.05)].sort_values("p_value")
    both     = sorted(set(kaks_pos["gene"]) & set(mk_sig["gene"]))
    n_all, n_ok = len(k), int((k["status"] == "OK").sum())
    print(f"Ka/Ks>1: {len(kaks_pos)} | MK sig: {len(mk_sig)} | convergents: {both}")

    # ----- 1. Distribution des Ka/Ks -----
    x = k["Ka_Ks"][(k["Ka_Ks"] > 0) & np.isfinite(k["Ka_Ks"])]
    fig, ax = plt.subplots(figsize=(7.2,4.4))
    bins = np.logspace(np.log10(x.min()), np.log10(x.max()), 46)
    ax.hist(x[x<=1], bins=bins, color=GREY, edgecolor="white", lw=.3, label="Ka/Ks ≤ 1")
    ax.hist(x[x>1],  bins=bins, color=RED,  edgecolor="white", lw=.3, label="Ka/Ks > 1 (candidats)")
    ax.axvline(1, color="black", ls="--", lw=1); ax.set_xscale("log")
    ax.set_xlabel("Ka/Ks (échelle logarithmique)"); ax.set_ylabel("Nombre de gènes")
    ax.set_title("Distribution des Ka/Ks sur les gènes à copie unique", color=NAVY, weight="bold")
    ax.legend(frameon=False); ax.spines[["top","right"]].set_visible(False)
    ax.text(.98,.95,f"{int((x>1).sum())} gènes avec Ka/Ks > 1", transform=ax.transAxes, ha="right", va="top", color=RED, weight="bold")
    fig.tight_layout(); fig.savefig(f"{out}/fig_distribution_kaks.png", dpi=300); plt.close(fig)

    # ----- 2. Ka vs Ks -----
    kk = k[(k["Ka"]>0) & (k["Ks"]>0) & (k["status"]=="OK")]
    ab = kk["Ka"] > kk["Ks"]
    fig, ax = plt.subplots(figsize=(6.6,6.2))
    ax.scatter(kk["Ks"][~ab], kk["Ka"][~ab], s=9, color=GREY, alpha=.5, lw=0, label="Ka/Ks ≤ 1")
    ax.scatter(kk["Ks"][ab],  kk["Ka"][ab],  s=12, color=RED, alpha=.7, lw=0, label="Ka/Ks > 1")
    lim=[min(kk["Ks"].min(),kk["Ka"].min()), max(kk["Ks"].max(),kk["Ka"].max())]
    ax.plot(lim, lim, "--", color="black", lw=1, label="Ka = Ks")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Ks (substitutions synonymes)"); ax.set_ylabel("Ka (substitutions non-synonymes)")
    ax.set_title("Ka en fonction de Ks", color=NAVY, weight="bold"); ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(); fig.savefig(f"{out}/fig_ka_vs_ks.png", dpi=300); plt.close(fig)

    # ----- 3. Volcano MK -----
    y = -np.log10(m["p_value"].clip(lower=1e-300)); sig = (m["p_value"]<0.05) & (m["alpha"]>0)
    fig, ax = plt.subplots(figsize=(7.2,4.7))
    ax.scatter(m["DoS"][~sig], y[~sig], s=14, color=GREY, alpha=.55, lw=0, label="non significatif")
    ax.scatter(m["DoS"][sig],  y[sig],  s=30, color=GREEN, edgecolor="black", lw=.3, label=f"α > 0 et p < 0,05 ({int(sig.sum())} gènes)", zorder=3)
    ax.axhline(-np.log10(0.05), color="black", ls="--", lw=1)
    ax.set_xlabel("DoS (Direction of Selection)"); ax.set_ylabel("-log10(p) — test de McDonald-Kreitman")
    ax.set_title("Test de McDonald-Kreitman (volcano plot)", color=NAVY, weight="bold")
    ax.legend(frameon=False, loc="upper left"); ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{out}/fig_volcano_mk.png", dpi=300); plt.close(fig)

    # ----- 4. Distribution du DoS -----
    dos = m["DoS"].dropna()
    fig, ax = plt.subplots(figsize=(7.2,4.2))
    ax.hist(dos, bins=60, color=GREY, edgecolor="white", lw=.2)
    ax.axvline(0, color="black", ls="--", lw=1)
    ax.set_xlabel("DoS (Direction of Selection)"); ax.set_ylabel("Nombre de gènes")
    ax.set_title("Distribution du DoS sur l'ensemble des gènes", color=NAVY, weight="bold")
    ax.spines[["top","right"]].set_visible(False)
    ax.text(.98,.95,f"médiane DoS = {dos.median():.3f}", transform=ax.transAxes, ha="right", va="top", color="#555")
    fig.tight_layout(); fig.savefig(f"{out}/fig_dos_distribution.png", dpi=300); plt.close(fig)

    # ----- 5. Entonnoir -----
    labels=["Gènes analysés","Ka/Ks calculable","Ka/Ks > 1","MK significatif","Convergents"]
    vals=[n_all,n_ok,len(kaks_pos),len(mk_sig),len(both)]; cols=[NAVY,BLUE,"#C0504D",GREEN,GOLD]
    fig, ax = plt.subplots(figsize=(7.6,4.4)); yy=np.arange(len(vals))[::-1]
    for yi,v,l,c in zip(yy,vals,labels,cols):
        w=np.log10(v+1); ax.barh(yi,w,color=c,height=0.62)
        ax.text(w+0.04,yi,f"{v}",va="center",weight="bold",color=c)
        ax.text(-0.05,yi,l,va="center",ha="right")
    ax.set_xlim(-0.05,5); ax.axis("off")
    ax.set_title("Entonnoir du pipeline", color=NAVY, weight="bold")
    fig.tight_layout(); fig.savefig(f"{out}/fig_funnel_pipeline.png", dpi=300); plt.close(fig)

    # ----- 6. Barplot 15 MK (par DoS) -----
    mm=mk_sig.sort_values("DoS")
    fig, ax = plt.subplots(figsize=(7.4,5.2))
    colors=[GOLD if g in both else GREEN for g in mm["gene"]]
    ax.barh(range(len(mm)), mm["DoS"], color=colors, edgecolor="white")
    ax.set_yticks(range(len(mm))); ax.set_yticklabels(mm["gene"], fontfamily="monospace", fontsize=9)
    for i,(d,p) in enumerate(zip(mm["DoS"],mm["p_value"])): ax.text(d+0.005,i,f"p={p:.3f}",va="center",fontsize=7.5,color="#555")
    ax.set_xlabel("DoS"); ax.set_title("Les gènes significatifs au test de McDonald-Kreitman", color=NAVY, weight="bold")
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{out}/fig_barplot_mk.png", dpi=300); plt.close(fig)

    # ----- tableaux (CSV) -----
    mk_sig[["gene","Dn","Ds","Pn","Ps","alpha","DoS","p_value","p_value_adjusted"]].to_csv(f"{out}/table_mk_candidats.csv", index=False)
    kaks_pos.sort_values("p_fisher").head(15)[["gene","Ka","Ks","Ka_Ks","p_fisher","p_fisher_adj"]].to_csv(f"{out}/table_kaks_top.csv", index=False)
    print("Terminé : 6 figures + 2 tableaux CSV écrits dans", out)

if __name__ == "__main__":
    a = sys.argv
    if len(a) < 3:
        print(__doc__)
    else:
        main(a[1], a[2], a[3] if len(a) > 3 else ".")
