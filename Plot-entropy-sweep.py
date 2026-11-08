import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from syndrome_entropy import measure


def sweep(distances, p_values, shots, max_bits, seed):
    results = {}
    for d in distances:
        rows = []
        for p in p_values:
            r = measure(d, float(p), shots, max_bits, seed)
            rows.append((p, r["H_faults"], r["H_marg"], r["H_cond"], r["unique_frac"]))
            print(f"d={d} p={p:.2e} H_faults={r['H_faults']:.4f} "
                  f"H_marg={r['H_marg']:.4f} H_cond={r['H_cond']:.4f} "
                  f"unique={r['unique_frac']:.3f}", flush=True)
        results[d] = np.array(rows)
    return results


def make_figure(results, out_path):
    fig, (ax_h, ax_r) = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, (d, arr) in enumerate(sorted(results.items())):
        c = colors[i % len(colors)]
        p, h_faults, h_marg, h_cond, uniq = arr.T
        ax_h.loglog(p, h_faults, "--", color=c, label=f"d={d}  H(P) bound")
        ax_h.loglog(p, h_cond, "o-", color=c, label=f"d={d}  H(s) measured")
        ax_r.semilogx(p, h_cond / h_faults, "o-", color=c, label=f"d={d}")
        ok = ~np.isnan(h_cond)
        if not ok.all():
            ax_h.plot([], [], " ", label=f"d={d}: p > {p[ok].max():.1e} dropped (undersampled)")

    ax_h.set_xlabel("physical error rate p")
    ax_h.set_ylabel("entropy per QEC round (nats)")
    ax_h.set_title("Syndrome entropy vs independent-fault bound")
    ax_h.grid(True, which="both", alpha=0.3)
    ax_h.legend(fontsize=8)

    ax_r.axhline(1.0, color="gray", linestyle=":", linewidth=1)
    ax_r.set_xlabel("physical error rate p")
    ax_r.set_ylabel("H(s) measured / H(P)")
    ax_r.set_title("How much H(s) ~ H(P) overstates the floor")
    ax_r.set_ylim(0, 1.1)
    ax_r.grid(True, which="both", alpha=0.3)
    ax_r.legend(fontsize=8)

    fig.savefig(out_path, dpi=200)
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    ap.add_argument("--p-min", type=float, default=3e-4)
    ap.add_argument("--p-max", type=float, default=1e-2)
    ap.add_argument("--points", type=int, default=6)
    ap.add_argument("--shots", type=int, default=2_000_000)
    ap.add_argument("--max-bits", type=int, default=64)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", default="entropy_sweep.png")
    args = ap.parse_args()

    p_values = np.logspace(np.log10(args.p_min), np.log10(args.p_max), args.points)
    results = sweep(args.distances, p_values, args.shots, args.max_bits, args.seed)
    make_figure(results, args.out)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
