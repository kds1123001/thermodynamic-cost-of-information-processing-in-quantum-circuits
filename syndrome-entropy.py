
import argparse
import math
from collections import Counter

import numpy as np
import stim

K_B = 1.380649e-23
HBAR = 1.054571817e-34


def h_nats(p: float) -> float:
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log(p) - (1.0 - p) * math.log(1.0 - p)


def build_circuit(basis: str, distance: int, rounds: int, p: float) -> stim.Circuit:
    return stim.Circuit.generated(
        f"surface_code:rotated_memory_{basis}",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )


def fault_entropy(circuit: stim.Circuit) -> float:
    total = 0.0

    def walk(c: stim.Circuit, mult: int) -> None:
        nonlocal total
        for inst in c:
            if isinstance(inst, stim.CircuitRepeatBlock):
                walk(inst.body_copy(), mult * inst.repeat_count)
                continue
            name = inst.name
            args = inst.gate_args_copy()
            targets = inst.targets_copy()
            if name == "DEPOLARIZE1":
                p = args[0]
                total += mult * len(targets) * (h_nats(p) + p * math.log(3))
            elif name == "DEPOLARIZE2":
                p = args[0]
                total += mult * (len(targets) // 2) * (h_nats(p) + p * math.log(15))
            elif name in ("X_ERROR", "Z_ERROR", "Y_ERROR"):
                p = args[0]
                total += mult * len(targets) * h_nats(p)
            elif name in ("MR", "MX", "MY", "M", "MRX", "MRY", "MRZ", "MZ"):

                if args and args[0] > 0:
                    total += mult * len(targets) * h_nats(args[0])
            elif name in ("R", "RX", "RY", "RZ"):
                pass
    walk(circuit, 1)
    return total


def per_round_fault_entropy(basis: str, distance: int, p: float) -> float:
    r = max(distance, 3)
    return fault_entropy(build_circuit(basis, distance, r + 1, p)) - fault_entropy(
        build_circuit(basis, distance, r, p)
    )


def plugin_entropy_mm(rows: np.ndarray) -> float:
    shots = rows.shape[0]
    packed = np.packbits(rows, axis=1)
    counts = Counter(map(bytes, packed))
    probs = np.array(list(counts.values()), dtype=np.float64) / shots
    h = float(-(probs * np.log(probs)).sum())
    return h + (len(counts) - 1) / (2.0 * shots)


def marginal_entropy_sum(rows: np.ndarray) -> float:
    freq = rows.mean(axis=0)
    return float(sum(h_nats(float(f)) for f in freq))


def measure(distance: int, p: float, shots: int, max_bits: int, seed: int,
            max_unique_frac: float = 0.1):
    rounds = max(distance, 3) + 3
    circuit = build_circuit("z", distance, rounds, p)
    by_t = {}
    for idx, coords in circuit.get_detector_coordinates().items():
        by_t.setdefault(coords[2], []).append(idx)
    times = sorted(by_t)
    t_mid = len(times) // 2
    cur = by_t[times[t_mid]]
    prev = by_t[times[t_mid - 1]]

    det = circuit.compile_detector_sampler(seed=seed).sample(shots=shots)
    rows_cur = det[:, cur]
    rows_pair = det[:, prev + cur]

    h_marg = marginal_entropy_sum(rows_cur)
    unique_frac = len(np.unique(np.packbits(rows_pair, axis=1), axis=0)) / shots
    if len(prev) + len(cur) <= max_bits and unique_frac <= max_unique_frac:
        h_cond = plugin_entropy_mm(rows_pair) - plugin_entropy_mm(det[:, prev])
    else:
        h_cond = float("nan")
    return {
        "bits": len(cur),
        "H_marg": h_marg,
        "H_cond": h_cond,
        "H_faults": per_round_fault_entropy("z", distance, p),
        "unique_frac": unique_frac,
    }


def landauer_power(h_nats_per_round: float, t_cold: float, t_hot: float, round_time: float):
    q_cold = K_B * t_cold * h_nats_per_round / round_time
    q_wall = K_B * (t_hot - t_cold) * h_nats_per_round / round_time
    return q_cold, q_wall


def reset_correction(m_ancillas: int, hbar_omega_over_k: float, t_cold: float) -> float:
    eps_r = 1.0 / (1.0 + math.exp(hbar_omega_over_k / t_cold))
    return m_ancillas * h_nats(eps_r)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", type=int, nargs="+", default=[3, 5])
    ap.add_argument("--p", type=float, default=1e-3, help="physical error rate")
    ap.add_argument("--shots", type=int, default=2_000_000)
    ap.add_argument("--max-bits", type=int, default=64,
                    help="skip the conditional estimate above this many detector bits (two rounds)")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--round-time", type=float, default=1e-6, help="seconds")
    ap.add_argument("--t-cold", type=float, default=0.020, help="K")
    ap.add_argument("--t-hot", type=float, default=300.0, help="K")
    ap.add_argument("--qubit-freq-ghz", type=float, default=5.0)
    args = ap.parse_args()

    hbar_omega_over_k = HBAR * 2 * math.pi * args.qubit_freq_ghz * 1e9 / K_B

    print(f"p = {args.p:g}, shots = {args.shots:,}, round time = {args.round_time:g} s")
    print(f"{'d':>3} {'bits/rd':>8} {'H_faults':>10} {'H_marg':>10} {'H_cond':>10} {'cond/faults':>12}")

    for d in args.distances:
        r = measure(d, args.p, args.shots, args.max_bits, args.seed)
        ratio = r["H_cond"] / r["H_faults"]
        print(f"{d:>3} {r['bits']:>8} {r['H_faults']:>10.4f} {r['H_marg']:>10.4f} "
              f"{r['H_cond']:>10.4f} {ratio:>12.4f}")

        h_used = r["H_cond"] if not math.isnan(r["H_cond"]) else r["H_faults"]
        q_cold, q_wall = landauer_power(h_used, args.t_cold, args.t_hot, args.round_time)
        corr = reset_correction(d * d - 1, hbar_omega_over_k, args.t_cold)
        print(f"      per logical qubit: cold {q_cold:.3e} W, wall-plug {q_wall:.3e} W, "
              f"reset correction m*h(eps_r) = {corr:.3e} nats")

    print("\nH_faults is the independent-fault bound H(P); H_cond is measured.")
    print("cond/faults < 1 quantifies how much H(s) ~ H(P) overstates the floor.")


if __name__ == "__main__":
    main()
