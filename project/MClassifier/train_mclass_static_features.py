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


def build_features(Xn, orders, mode, structure=None):
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
        static = np.broadcast_to(structure.reshape(1, 1, -1), (B, M, structure.size))
        parts.append(static.astype(np.float32))

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
    parser.add_argument("--mode", choices=["prefix", "seq_prefix", "signed_seq_prefix", "signed_seq_prefix_structure"], default="prefix")
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

    idx = np.arange(len(argmin))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    Xn, y_mean, y_std = standardize_y(X[train_idx], X)
    F_all = build_features(Xn, orders, args.mode, structure=structure)
    F_train, F_val, f_mean, f_std = standardize_features(F_all[train_idx], F_all[val_idx])
    y_train = basin[train_idx].astype(np.float32)
    y_val = basin[val_idx].astype(np.float32)
    argmin_val = argmin[val_idx].astype(np.int64)

    input_dim = F_train.shape[-1]
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
    loader = DataLoader(
        TensorDataset(
            torch.from_numpy(F_train.reshape(-1, input_dim)),
            torch.from_numpy(y_train.reshape(-1)),
        ),
        batch_size=args.batch_rows,
        shuffle=True,
    )

    Fv = torch.from_numpy(F_val.astype(np.float32))
    yv = torch.from_numpy(y_val.astype(np.float32))
    argmin_v = torch.from_numpy(argmin_val)
    pos_w = torch.tensor(pos_weight, dtype=torch.float32)

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        for xb, yb in loader:
            opt.zero_grad()
            logits = model(xb).reshape(-1)
            loss = F.binary_cross_entropy_with_logits(logits, yb, pos_weight=pos_w)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * xb.shape[0]
            total += xb.shape[0]

        model.eval()
        with torch.no_grad():
            Bv = Fv.shape[0]
            scores = model(Fv.reshape(Bv * M, input_dim)).reshape(Bv, M)
            val_loss = float(F.binary_cross_entropy_with_logits(scores, yv, pos_weight=pos_w).item())
            ranked = torch.topk(scores, k=M, dim=1, largest=True).indices
            row = torch.arange(Bv)
            basin1 = float(yv[row, ranked[:, 0]].mean().item())
            basin4 = float(yv.gather(1, ranked[:, :min(4, M)]).amax(dim=1).mean().item())
            basin8 = float(yv.gather(1, ranked[:, :min(8, M)]).amax(dim=1).mean().item())
            arg1 = float((ranked[:, 0] == argmin_v).float().mean().item())
            arg8 = float((ranked[:, :min(8, M)] == argmin_v[:, None]).any(dim=1).float().mean().item())
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
