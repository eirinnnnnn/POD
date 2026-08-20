#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_basin import load_basin_dataset, standardize


def to_llr(X, meta):
    snr = float(meta.get("awgn_start", 2.0))
    rate = float(meta.get("k", 1)) / float(meta.get("n", X.shape[1]))
    var = (10.0 ** (-snr / 10.0)) / rate / 2.0
    return X * (2.0 / var)


def main():
    parser = argparse.ArgumentParser(description="Train basin classifier with LLR input.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--delta", type=float, default=1e-12)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
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
    X = to_llr(X, meta)

    idx = np.arange(len(argmin))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    if args.no_standardize:
        X_train = X[train_idx].astype(np.float32)
        X_val = X[val_idx].astype(np.float32)
        mean = np.zeros((1, X.shape[1]), dtype=np.float32)
        std = np.ones((1, X.shape[1]), dtype=np.float32)
    else:
        X_train, X_val, mean, std = standardize(X[train_idx], X[val_idx])

    y_train = basin[train_idx]
    y_val = basin[val_idx]
    argmin_val = argmin[val_idx]

    pos_rate = float(y_train.mean())
    pos_weight = (1.0 - pos_rate) / max(pos_rate, 1e-6)

    model = nn.Sequential(
        nn.Linear(X.shape[1], args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, args.hidden),
        nn.ReLU(),
        nn.Linear(args.hidden, M),
    )
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.full((M,), pos_weight))

    train_loader = DataLoader(
        TensorDataset(torch.from_numpy(X_train), torch.from_numpy(y_train.astype(np.float32))),
        batch_size=args.batch,
        shuffle=True,
    )

    Xv = torch.from_numpy(X_val.astype(np.float32))
    yv = torch.from_numpy(y_val.astype(np.float32))
    argmin_v = torch.from_numpy(argmin_val.astype(np.int64))

    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        for xb, yb in train_loader:
            opt.zero_grad()
            logits = model(xb)
            loss = loss_fn(logits, yb)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * len(yb)
            total += len(yb)

        model.eval()
        with torch.no_grad():
            logits = model(Xv)
            val_loss = float(loss_fn(logits, yv).item())
            probs = torch.sigmoid(logits)
            top1 = probs.argmax(dim=1)
            top4 = probs.topk(k=min(4, M), dim=1).indices
            top8 = probs.topk(k=min(8, M), dim=1).indices
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
            "kind": "multi_label_best_basin_classifier_llr",
            "state_dict": model.state_dict(),
            "n": X.shape[1],
            "M": M,
            "hidden": args.hidden,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "delta": args.delta,
            "input_mode": "llr",
            "standardize": (not args.no_standardize),
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
