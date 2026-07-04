#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_branch_ranker import load_orders, permute_batch
from train_mclass_ranker import load_metric_dataset, standardize


def target_scale(gaps):
    positive = gaps[gaps > 1e-9]
    if positive.size == 0:
        return 1.0
    scale = float(np.median(positive))
    return max(scale, 1e-9)


def make_targets(metrics, scale, mode):
    gaps = metrics - metrics.min(axis=1, keepdims=True)
    if mode == "log_gap":
        return np.log1p(gaps / scale).astype(np.float32)
    if mode == "gap":
        return (gaps / scale).astype(np.float32)
    raise ValueError(mode)


def main():
    parser = argparse.ArgumentParser(description="Train shared score(permuted-y) -> final SCL metric gap.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--orders", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=64)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--target", choices=["log_gap", "gap"], default="log_gap")
    parser.add_argument("--loss", choices=["mse", "huber"], default="huber")
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    rng = np.random.default_rng(args.seed)
    X, metrics, label, M, meta = load_metric_dataset(args.dataset)
    orders = load_orders(args.orders)
    if orders.shape != (M, X.shape[1]):
        raise SystemExit(f"orders shape {orders.shape} incompatible with n={X.shape[1]} M={M}")

    idx = np.arange(len(label))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    X_train, X_val, mean, std = standardize(X[train_idx], X[val_idx])
    Xn = np.empty_like(X, dtype=np.float32)
    Xn[train_idx] = X_train
    Xn[val_idx] = X_val

    gaps = metrics - metrics.min(axis=1, keepdims=True)
    scale = target_scale(gaps[train_idx])
    targets = make_targets(metrics, scale, args.target)

    P_train = permute_batch(Xn[train_idx], orders)
    P_val = permute_batch(Xn[val_idx], orders)
    t_train = targets[train_idx]
    t_val = targets[val_idx]
    y_val = label[val_idx]

    model = nn.Sequential(
        nn.Linear(X.shape[1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, 1),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    loader = DataLoader(
        TensorDataset(torch.from_numpy(P_train), torch.from_numpy(t_train.astype(np.float32))),
        batch_size=args.batch,
        shuffle=True,
    )

    Pv = torch.from_numpy(P_val.astype(np.float32))
    tv = torch.from_numpy(t_val.astype(np.float32))
    yv = torch.from_numpy(y_val.astype(np.int64))

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        for xb, tb in loader:
            opt.zero_grad()
            B = xb.shape[0]
            pred = model(xb.reshape(B * M, -1)).reshape(B, M)
            if args.loss == "mse":
                loss = F.mse_loss(pred, tb)
            else:
                loss = F.smooth_l1_loss(pred, tb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * B
            total += B

        model.eval()
        with torch.no_grad():
            Bv = Pv.shape[0]
            val_pred = model(Pv.reshape(Bv * M, -1)).reshape(Bv, M)
            if args.loss == "mse":
                val_loss = float(F.mse_loss(val_pred, tv).item())
            else:
                val_loss = float(F.smooth_l1_loss(val_pred, tv).item())
            pred_argmin = val_pred.argmin(dim=1)
            top4 = val_pred.topk(k=min(4, M), dim=1, largest=False).indices
            top8 = val_pred.topk(k=min(8, M), dim=1, largest=False).indices
            acc1 = float((pred_argmin == yv).float().mean().item())
            acc4 = float((top4 == yv[:, None]).any(dim=1).float().mean().item())
            acc8 = float((top8 == yv[:, None]).any(dim=1).float().mean().item())

        print(
            f"epoch={epoch} train_{args.loss}={total_loss/total:.6f} val_{args.loss}={val_loss:.6f} "
            f"val_argmin_top1={acc1:.6f} val_argmin_top4={acc4:.6f} val_argmin_top8={acc8:.6f}"
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    import torch
    torch.save(
        {
            "kind": "shared_permuted_metric_gap_regressor",
            "state_dict": model.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "orders": orders,
            "target": args.target,
            "loss": args.loss,
            "target_scale": scale,
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
