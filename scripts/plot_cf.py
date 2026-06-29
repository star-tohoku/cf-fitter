#!/usr/bin/env python3
"""Plot phi-p and phi-alpha correlation functions (LL solid, KP dashed)."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(fname):
    with open(fname) as f:
        rows = list(csv.reader(f))
    head = rows[0]
    data = {h: [float(r[i]) for r in rows[1:]] for i, h in enumerate(head)}
    return data

# ---------------- phi-p ----------------
d = load("cf_phip.csv")
k = d["k_MeV"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
sets = [("HAL",          "tab:blue",  r"$^4S_{3/2}$ HAL QCD-like ($f_0$=1.43, $d_0$=2.35 fm)"),
        ("HALplusBound", "tab:red",   r"+ bound $^2S_{1/2}$ ($f_0$=$-$4.24 fm, B.E.$\sim$3 MeV)"),
        ("ALICEeff",     "tab:green", r"eff. complex $f_0$=0.85+0.16$i$ fm (LL only)")]
for ax, R in zip(axes, [1.2, 3.0]):
    for name, col, lab in sets:
        ax.plot(k, d[f"LL_{name}_R{R}"], color=col, lw=1.8, label=lab)
        key = f"KP_{name}_R{R}"
        if key in d:
            ax.plot(k, d[key], color=col, lw=1.8, ls="--")
    ax.axhline(1.0, color="gray", lw=0.6)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_title(rf"$\phi\,p$,  $R$ = {R} fm")
    ax.set_xlim(0, 300)
axes[0].set_ylabel(r"$C(k^*)$")
axes[0].legend(fontsize=8, title="solid: LL,  dashed: KP", loc="upper right")
fig.tight_layout()
fig.savefig("cf_phip.png", dpi=160)
fig.savefig("cf_phip.pdf")

# ---------------- phi-alpha ----------------
d = load("cf_phialpha.csv")
k = d["k_MeV"]
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
sets = [("weak",   "tab:blue",   r"weak: $V_0$=$-$6.5 MeV ($f_0$=1.54 fm)"),
        ("strong", "tab:orange", r"near-unitary: $V_0$=$-$15.0 MeV ($f_0$=24.5 fm)"),
        ("bound",  "tab:red",    r"bound: $V_0$=$-$21.9 MeV ($f_0$=$-$8.5 fm, B.E.$\sim$0.5 MeV)")]
for ax, R in zip(axes, [1.2, 2.5, 5.0]):
    for name, col, lab in sets:
        ax.plot(k, d[f"LL_{name}_R{R}"], color=col, lw=1.8, label=lab)
        ax.plot(k, d[f"KP_{name}_R{R}"], color=col, lw=1.8, ls="--")
    ax.axhline(1.0, color="gray", lw=0.6)
    ax.set_xlabel(r"$k^*$ [MeV/$c$]")
    ax.set_title(rf"$\phi\,\alpha$,  $R$ = {R} fm")
    ax.set_xlim(0, 200)
axes[0].set_ylabel(r"$C(k^*)$")
axes[0].set_yscale("log")
axes[0].legend(fontsize=8, title="solid: LL,  dashed: KP", loc="upper right")
fig.tight_layout()
fig.savefig("cf_phialpha.png", dpi=160)
fig.savefig("cf_phialpha.pdf")
print("plots written")
