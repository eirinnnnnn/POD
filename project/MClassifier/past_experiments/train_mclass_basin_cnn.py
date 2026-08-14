#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_basin import load_basin_dataset, standardize


def main():
    parser = argparse.ArgumentParser(description="Train 1D-CNN multi-label best-basin M-classifier.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--delta", type=float, default=1e-12)
    parser.add_argument("--channels", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=25)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--batch", type=int, default=256)
    parser.add_argument("--val-frac", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset
    except ImportError as exc:
        raise SystemExit("PyTorch is required for training: pip install torch") from exc

    class BasinCNN(nn.Module):
        def __init__(self, n, M, channels):
            super().__init__()
            self.conv = nn.Sequential(
                nn.Conv1d(1, channels, kernel_size=5, padding=2),
                nn.BatchNorm1d(channels),
                nn.ReLU(),
                nn.Conv1d(channels, channels, kernel_size=5, padding=2),
                nn.BatchNorm1d(channels),
                nn.ReLU(),
                nn.MaxPool1d(2),
                nn.Conv1d(channels, channels * 2, kernel_size=5, padding=2),
                nn.BatchNorm1d(channels * 2),
                nn.ReLU(),
                nn.Conv1d(channels * 2, channels * 2, kernel_size=3, padding=1),
                nn.BatchNorm1d(channels * 2),
                nn.ReLU(),
                nn.MaxPool1d(2),
            )
            self.head = nn.Sequential(
                nn.Flatten(),
                nn.Linear((channels * 2) * (n // 4), 256),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.Linear(256, M),
            )

        def forward(self, x):
            return self.head(self.conv(x[:, None, :]))

    rng = np.random.default_rng(args.seed)
    X, basin, argmin, _, M, meta = load_basin_dataset(args.dataset, args.delta)
    idx = np.arange(len(argmin))
    rng.shuffle(idx)
    val_count = max(1, int(len(idx) * args.val_frac))
    val_idx = idx[:val_count]
    train_idx = idx[val_count:]

    X_train, X_val, mean, std = standardize(X[train_idx], X[val_idx])
    y_train = basin[train_idx]
    y_val = basin[val_idx]
    argmin_val = argmin[val_idx]

    pos_rate = float(y_train.mean())
    pos_weight = (1.0 - pos_rate) / max(pos_rate, 1e-6)

    model = BasinCNN(X.shape[1], M, args.channels)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-5)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.full((M,), pos_weight))
    loader = DataLoader(
        TensorDataset(torch.from_numpy(X_train.astype(np.float32)), torch.from_numpy(y_train.astype(np.float32))),
        batch_size=args.batch,
        shuffle=True,
    )

    Xv = torch.from_numpy(X_val.astype(np.float32))
    yv = torch.from_numpy(y_val.astype(np.float32))
    argmin_v = torch.from_numpy(argmin_val.astype(np.int64))

    best_state = None
    best_score = -1.0
    best_epoch = 0
    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        total = 0
        for xb, yb in loader:
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
            top1 = logits.argmax(dim=1)
            top4 = logits.topk(k=min(4, M), dim=1).indices
            top8 = logits.topk(k=min(8, M), dim=1).indices
            row = torch.arange(len(yv))
            basin_top1 = float(yv[row, top1].mean().item())
            basin_top4 = float(yv.gather(1, top4).max(dim=1).values.mean().item())
            basin_top8 = float(yv.gather(1, top8).max(dim=1).values.mean().item())
            argmin_top8 = float((top8 == argmin_v[:, None]).any(dim=1).float().mean().item())

        if basin_top8 > best_score:
            best_score = basin_top8
            best_epoch = epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}

        print(
            f"epoch={epoch} train_bce={total_loss/total:.6f} val_bce={val_loss:.6f} "
            f"basin_top1={basin_top1:.6f} basin_top4={basin_top4:.6f} "
            f"basin_top8={basin_top8:.6f} argmin_top8={argmin_top8:.6f}"
        )

    if best_state is not None:
        model.load_state_dict(best_state)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "kind": "cnn_multi_label_best_basin_classifier",
            "state_dict": model.state_dict(),
            "n": X.shape[1],
            "M": M,
            "channels": args.channels,
            "mean": mean.astype(np.float32),
            "std": std.astype(np.float32),
            "delta": args.delta,
            "best_epoch": best_epoch,
            "best_basin_top8": best_score,
            "meta": meta,
        },
        out,
    )
    print(f"wrote {out} best_epoch={best_epoch} best_basin_top8={best_score:.6f}")


if __name__ == "__main__":
    main()
