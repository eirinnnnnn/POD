#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_basin import load_basin_dataset


def main():
    parser = argparse.ArgumentParser(description="Predict top-k branches from CNN basin classifier.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--topk", type=int, default=1)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
    except ImportError as exc:
        raise SystemExit("PyTorch is required for prediction: pip install torch") from exc

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

    ckpt = torch.load(args.model, map_location="cpu")
    X, _, _, _, _, _ = load_basin_dataset(args.dataset, float(ckpt.get("delta", 1e-12)))
    model = BasinCNN(int(ckpt["n"]), int(ckpt["M"]), int(ckpt["channels"]))
    model.load_state_dict(ckpt["state_dict"])
    model.eval()

    Xn = ((X - ckpt["mean"]) / ckpt["std"]).astype(np.float32)
    with torch.no_grad():
        logits = model(torch.from_numpy(Xn))
        topk = torch.topk(logits, k=min(args.topk, int(ckpt["M"])), dim=1).indices.cpu().numpy()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("sample")
        for i in range(topk.shape[1]):
            f.write(f",idx{i}")
        f.write("\n")
        for s, row in enumerate(topk):
            f.write(str(s))
            for idx in row:
                f.write(f",{int(idx)}")
            f.write("\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
