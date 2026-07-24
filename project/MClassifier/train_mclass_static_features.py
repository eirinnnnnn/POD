#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_basin import load_basin_dataset
from train_mclass_branch_ranker import load_orders, permute_batch


PREFIXES = (8, 16, 32, 64, 128)


def load_structure(path, n):
    if not path:
        return None
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=np.float32)
    if data.shape == ():
        data = data.reshape(1)
    bit_type = data["bit_type"].astype(np.int64)
    rel_size = data["relation_size"].astype(np.float32)
    if len(bit_type) != n:
        raise ValueError(f"structure length {len(bit_type)} != n={n}")
    type_onehot = np.eye(3, dtype=np.float32)[np.clip(bit_type, 0, 2)]
    rel_norm = rel_size / max(float(rel_size.max()), 1.0)
    return np.concatenate([type_onehot, rel_norm[:, None]], axis=1).astype(np.float32)


def load_gmatrix(path, n):
    """Load the full generator matrix, transposed to (n, k): row i is the
    length-k support vector of message bits feeding codeword position i
    (as written by mclass_gmatrix_dump). This is the actual G, not a
    scalar summary of it."""
    if not path:
        return None
    raw = np.genfromtxt(path, delimiter=",", skip_header=1, dtype=np.float32)
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    gmat = raw[:, 1:]  # drop the leading "pos" column
    if gmat.shape[0] != n:
        raise ValueError(f"gmatrix length {gmat.shape[0]} != n={n}")
    return gmat.astype(np.float32)


def prefix_features(abs_perm, prefixes, thresholds):
    B, M, n = abs_perm.shape
    feats = []
    for p in prefixes:
        p = min(p, n)
        x = abs_perm[:, :, :p]
        feats.extend([
            x.min(axis=2),
            x.mean(axis=2),
            x.max(axis=2),
            x.std(axis=2),
        ])
        k = min(4, p)
        feats.append(np.partition(x, kth=k - 1, axis=2)[:, :, :k].mean(axis=2))
        for tau in thresholds:
            feats.append((x < tau).mean(axis=2))
    return np.stack(feats, axis=2).astype(np.float32)


def build_features(Xn, orders, mode, structure=None, gmatrix=None):
    perm = permute_batch(Xn, orders)
    abs_perm = np.abs(perm).astype(np.float32)
    sign_perm = np.sign(perm).astype(np.float32)
    prefixes = tuple(p for p in PREFIXES if p <= Xn.shape[1])
    if Xn.shape[1] not in prefixes:
        prefixes = prefixes + (Xn.shape[1],)

    parts = []
    if "seq" in mode:
        parts.append(abs_perm)
    if "signed" in mode:
        parts.append(perm)
        parts.append(sign_perm)
    if "prefix" in mode:
        parts.append(prefix_features(abs_perm, prefixes, thresholds=(0.25, 0.5, 1.0)))
    if "structure" in mode:
        if structure is None:
            raise ValueError("mode requested structure but no --structure was provided")
        B, M, _ = abs_perm.shape
        static_by_branch = structure[orders].reshape(1, M, -1)
        static = np.broadcast_to(static_by_branch, (B, M, static_by_branch.shape[-1]))
        parts.append(static.astype(np.float32))
    if "gmatrix" in mode:
        if gmatrix is None:
            raise ValueError("mode requested gmatrix but no --gmatrix was provided")
        B, M, _ = abs_perm.shape
        # gmatrix[i] = support vector (over k message bits) of codeword position i.
        # Permuting rows by `orders` mirrors exactly how the received vector
        # itself is permuted per branch, so this is literally the branch's
        # decoding-order view of G -- known before any decoding happens.
        gmat_by_branch = gmatrix[orders].reshape(1, M, -1)
        gmat = np.broadcast_to(gmat_by_branch, (B, M, gmat_by_branch.shape[-1]))
        parts.append(gmat.astype(np.float32))

    if not parts:
        raise ValueError(f"empty feature mode: {mode}")
    return np.concatenate([p.reshape(p.shape[0], p.shape[1], -1) for p in parts], axis=2)


def standardize_y(train_X, all_X):
    mean = train_X.mean(axis=0, keepdims=True)
    std = train_X.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return ((all_X - mean) / std).astype(np.float32), mean.astype(np.float32), std.astype(np.float32)


def standardize_features(train_F, val_F):
    flat = train_F.reshape(-1, train_F.shape[-1])
    mean = flat.mean(axis=0, keepdims=True)
    std = flat.std(axis=0, keepdims=True)
    std[std < 1e-6] = 1.0
    return ((train_F - mean) / std).astype(np.float32), ((val_F - mean) / std).astype(np.float32), mean.reshape(-1), std.reshape(-1)


def main():
    parser = argparse.ArgumentParser(description="Train hardware-realistic static branch features -> basin scorer.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--orders", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--structure", default="")
    parser.add_argument("--gmatrix", default="")
    parser.add_argument("--mode", choices=["prefix", "seq_prefix", "signed_seq_prefix", "signed_seq_prefix_structure", "signed_seq_prefix_gmatrix"], default="prefix")
    parser.add_argument("--delta", type=float, default=1e-12)
    parser.add_argument("--hidden", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch-rows", type=int, default=4096)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--max-samples", type=int, default=0)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, basin, argmin, _, M, meta = load_basin_dataset(args.dataset, args.delta)
    if args.max_samples and args.max_samples < len(argmin):
        X = X[:args.max_samples]
        basin = basin[:args.max_samples]
        argmin = argmin[:args.max_samples]
    orders = load_orders(args.orders)
    if orders.shape != (M, X.shape[1]):
        raise SystemExit(f"orders shape {orders.shape} incompatible with n={X.shape[1]} M={M}")
    structure = load_structure(args.structure, X.shape[1])
    gmatrix = load_gmatrix(args.gmatrix, X.shape[1])

    idx = np.arange(len(argmin))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    # "gmatrix" mode's static per-branch feature (M, n*k) is identical across
    # every sample, so materializing a (samples, M, n*k) broadcast copy is
    # wasteful and can blow past available RAM (e.g. 6000x64x8192 floats =
    # 11.7 GiB). Instead keep the sample-varying (LLR-derived) features and
    # the branch-only static features separate, and only concatenate them
    # for the rows actually needed by the current mini-batch.
    heavy_static = gmatrix is not None and "gmatrix" in args.mode
    dyn_mode = args.mode.replace("_gmatrix", "") if heavy_static else args.mode

    Xn, y_mean, y_std = standardize_y(X[train_idx], X)
    F_dyn_all = build_features(Xn, orders, dyn_mode, structure=structure)
    F_dyn_train, F_dyn_val, f_mean, f_std = standardize_features(F_dyn_all[train_idx], F_dyn_all[val_idx])
    y_train = basin[train_idx].astype(np.float32)
    y_val = basin[val_idx].astype(np.float32)
    argmin_val = argmin[val_idx].astype(np.int64)

    dyn_dim = F_dyn_train.shape[-1]
    static_g = gmatrix[orders].reshape(M, -1).astype(np.float32) if heavy_static else None
    static_dim = static_g.shape[-1] if heavy_static else 0
    input_dim = dyn_dim + static_dim

    model = nn.Sequential(
        nn.Linear(input_dim, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    pos_rate = float(y_train.mean())
    pos_weight = (1.0 - pos_rate) / max(pos_rate, 1e-6)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)
    pos_w = torch.tensor(pos_weight, dtype=torch.float32)
    static_g_t = torch.from_numpy(static_g) if heavy_static else None

    def make_batch(dyn_np, gt_np=None):
        xb = torch.from_numpy(dyn_np)
        b = xb.shape[0]
        if static_g_t is not None:
            xb = torch.cat([xb, static_g_t.unsqueeze(0).expand(b, M, static_dim)], dim=2)
        xb = xb.reshape(b * M, input_dim)
        if gt_np is None:
            return xb
        yb = torch.from_numpy(gt_np).reshape(b * M)
        return xb, yb

    samples_per_batch = max(1, args.batch_rows // M)

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        perm_idx = rng.permutation(F_dyn_train.shape[0])
        for start in range(0, len(perm_idx), samples_per_batch):
            sel = perm_idx[start:start + samples_per_batch]
            xb, yb = make_batch(F_dyn_train[sel], y_train[sel])
            opt.zero_grad()
            logits = model(xb).reshape(-1)
            loss = F.binary_cross_entropy_with_logits(logits, yb, pos_weight=pos_w)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * xb.shape[0]
            total += xb.shape[0]

        model.eval()
        with torch.no_grad():
            Bv = F_dyn_val.shape[0]
            val_loss_sum = 0.0
            basin1_hits = []
            basin4_hits = []
            basin8_hits = []
            arg1_hits = []
            arg8_hits = []
            for start in range(0, Bv, samples_per_batch):
                sel = np.arange(start, min(start + samples_per_batch, Bv))
                xb, yb = make_batch(F_dyn_val[sel], y_val[sel])
                b = len(sel)
                scores = model(xb).reshape(b, M)
                yv_chunk = yb.reshape(b, M)
                val_loss_sum += float(F.binary_cross_entropy_with_logits(scores, yv_chunk, pos_weight=pos_w).item()) * b
                ranked = torch.topk(scores, k=M, dim=1, largest=True).indices
                row = torch.arange(b)
                argmin_chunk = torch.from_numpy(argmin_val[sel])
                basin1_hits.append(yv_chunk[row, ranked[:, 0]])
                basin4_hits.append(yv_chunk.gather(1, ranked[:, :min(4, M)]).amax(dim=1))
                basin8_hits.append(yv_chunk.gather(1, ranked[:, :min(8, M)]).amax(dim=1))
                arg1_hits.append((ranked[:, 0] == argmin_chunk).float())
                arg8_hits.append((ranked[:, :min(8, M)] == argmin_chunk[:, None]).any(dim=1).float())
            val_loss = val_loss_sum / Bv
            basin1 = float(torch.cat(basin1_hits).mean().item())
            basin4 = float(torch.cat(basin4_hits).mean().item())
            basin8 = float(torch.cat(basin8_hits).mean().item())
            arg1 = float(torch.cat(arg1_hits).mean().item())
            arg8 = float(torch.cat(arg8_hits).mean().item())
        print(
            f"epoch={epoch} train_bce={total_loss/total:.6f} val_bce={val_loss:.6f} "
            f"basin_top1={basin1:.6f} basin_top4={basin4:.6f} basin_top8={basin8:.6f} "
            f"argmin_top1={arg1:.6f} argmin_top8={arg8:.6f}",
            flush=True,
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    import torch
    torch.save(
        {
            "kind": "static_branch_feature_basin_scorer",
            "state_dict": model.state_dict(),
            "mode": args.mode,
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "input_dim": input_dim,
            "orders": orders,
            "structure": structure,
            "gmatrix": gmatrix,
            "y_mean": y_mean,
            "y_std": y_std,
            "feature_mean": f_mean.astype(np.float32),
            "feature_std": f_std.astype(np.float32),
            "delta": args.delta,
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
