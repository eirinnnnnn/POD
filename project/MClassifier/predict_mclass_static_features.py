#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass_basin import load_basin_dataset
from train_mclass_static_features import build_features


def main():
    parser = argparse.ArgumentParser(description="Predict top-k branches from static branch feature scorer.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--topk", type=int, default=64)
    parser.add_argument("--batch-samples", type=int, default=256)
    args = parser.parse_args()

    try:
        import torch
        import torch.nn as nn
    except ImportError as exc:
        raise SystemExit("PyTorch is required for prediction: pip install torch") from exc

    ckpt = torch.load(args.model, map_location="cpu")
    X, _, _, _, M, _ = load_basin_dataset(args.dataset, float(ckpt.get("delta", 1e-12)))
    if M != int(ckpt["M"]):
        raise SystemExit(f"dataset M={M} but model M={ckpt['M']}")

    Xn = ((X - ckpt["y_mean"]) / ckpt["y_std"]).astype(np.float32)
    features = build_features(
        Xn,
        np.asarray(ckpt["orders"], dtype=np.int64),
        str(ckpt["mode"]),
        structure=ckpt.get("structure"),
    )
    mean = np.asarray(ckpt["feature_mean"], dtype=np.float32).reshape(1, 1, -1)
    std = np.asarray(ckpt["feature_std"], dtype=np.float32).reshape(1, 1, -1)
    features = ((features - mean) / std).astype(np.float32)

    model = nn.Sequential(
        nn.Linear(int(ckpt["input_dim"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), 1),
    )
    model.load_state_dict(ckpt["state_dict"])
    model.eval()

    top_rows = []
    with torch.no_grad():
        for start in range(0, features.shape[0], args.batch_samples):
            xb = torch.from_numpy(features[start:start + args.batch_samples])
            B = xb.shape[0]
            scores = model(xb.reshape(B * M, -1)).reshape(B, M)
            top = torch.topk(scores, k=min(args.topk, M), dim=1, largest=True).indices.cpu().numpy()
            top_rows.append(top)
    topk = np.concatenate(top_rows, axis=0)

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
