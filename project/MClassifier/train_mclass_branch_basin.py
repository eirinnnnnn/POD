#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_basin import load_basin_dataset, standardize
from train_mclass_basin_llr import to_llr
from train_mclass_branch_ranker import load_orders, permute_batch


def main():
    parser = argparse.ArgumentParser(description="Train shared score(permuted-y) best-basin classifier.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--orders", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--delta", type=float, default=1e-12)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=3)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--input-mode", choices=["raw", "llr"], default="raw")
    parser.add_argument("--no-standardize", action="store_true")
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, basin, argmin, _, M, meta = load_basin_dataset(args.dataset, args.delta)
    orders = load_orders(args.orders)
    if orders.shape[0] != M or orders.shape[1] != X.shape[1]:
        raise SystemExit(f"orders shape {orders.shape} incompatible with dataset n={X.shape[1]} M={M}")
    if args.input_mode == "llr":
        X = to_llr(X, meta)

    idx = np.arange(len(argmin))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    if args.no_standardize:
        Xn = X.astype(np.float32)
        mean = np.zeros((1, X.shape[1]), dtype=np.float32)
        std = np.ones((1, X.shape[1]), dtype=np.float32)
    else:
        X_train, X_val, mean, std = standardize(X[train_idx], X[val_idx])
        Xn = np.empty_like(X, dtype=np.float32)
        Xn[train_idx] = X_train
        Xn[val_idx] = X_val

    P_train = permute_batch(Xn[train_idx], orders)
    P_val = permute_batch(Xn[val_idx], orders)
    y_train = basin[train_idx].astype(np.float32)
    y_val = basin[val_idx].astype(np.float32)
    argmin_val = argmin[val_idx].astype(np.int64)

    pos_rate = float(y_train.mean())
    pos_weight = (1.0 - pos_rate) / max(pos_rate, 1e-6)

    scorer = nn.Sequential(
        nn.Linear(X.shape[1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    opt = torch.optim.Adam(scorer.parameters(), lr=args.lr, weight_decay=1e-5)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.full((M,), pos_weight))
    loader = DataLoader(
        TensorDataset(torch.from_numpy(P_train), torch.from_numpy(y_train)),
        batch_size=args.batch,
        shuffle=True,
    )

    Pv = torch.from_numpy(P_val)
    yv = torch.from_numpy(y_val)
    argmin_v = torch.from_numpy(argmin_val)

    for epoch in range(1, args.epochs + 1):
        scorer.train()
        total_loss = 0.0
        total = 0
        for xb, yb in loader:
            opt.zero_grad()
            B = xb.shape[0]
            logits = scorer(xb.reshape(B * M, -1)).reshape(B, M)
            loss = loss_fn(logits, yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * B
            total += B

        scorer.eval()
        with torch.no_grad():
            Bv = Pv.shape[0]
            val_logits = scorer(Pv.reshape(Bv * M, -1)).reshape(Bv, M)
            val_loss = float(loss_fn(val_logits, yv).item())
            top1 = val_logits.argmax(dim=1)
            top4 = val_logits.topk(k=min(4, M), dim=1).indices
            top8 = val_logits.topk(k=min(8, M), dim=1).indices
            row = torch.arange(len(yv))
            basin_top1 = float(yv[row, top1].mean().item())
            basin_top4 = float(yv.gather(1, top4).max(dim=1).values.mean().item())
            basin_top8 = float(yv.gather(1, top8).max(dim=1).values.mean().item())
            argmin_top1 = float((top1 == argmin_v).float().mean().item())
            argmin_top8 = float((top8 == argmin_v[:, None]).any(dim=1).float().mean().item())

        print(
            f"epoch={epoch} train_bce={total_loss/total:.6f} val_bce={val_loss:.6f} "
            f"basin_top1={basin_top1:.6f} basin_top4={basin_top4:.6f} basin_top8={basin_top8:.6f} "
            f"argmin_top1={argmin_top1:.6f} argmin_top8={argmin_top8:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "kind": "shared_permuted_branch_basin_classifier",
            "state_dict": scorer.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "orders": orders,
            "delta": args.delta,
            "input_mode": args.input_mode,
            "standardize": (not args.no_standardize),
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
