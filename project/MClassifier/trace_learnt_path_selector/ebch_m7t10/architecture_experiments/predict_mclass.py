#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np

from train_mclass import load_dataset


def main():
    parser = argparse.ArgumentParser(description="Predict top-k M-classifier branches.")
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

    X, _, _, _ = load_dataset(args.dataset)
    ckpt = torch.load(args.model, map_location="cpu")
    model = nn.Sequential(
        nn.Linear(int(ckpt["n"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), int(ckpt["hidden"])),
        nn.ReLU(),
        nn.Linear(int(ckpt["hidden"]), int(ckpt["M"])),
    )
    model.load_state_dict(ckpt["state_dict"])
    model.eval()

    Xn = (X - ckpt["mean"]) / ckpt["std"]
    with torch.no_grad():
        logits = model(torch.from_numpy(Xn.astype(np.float32)))
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
