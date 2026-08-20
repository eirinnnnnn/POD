# Adaptive Path Selector -- plan and methodology

Code: m7t10_m64 (eBCH, n=128, k=64), AED M=64 automorphisms, SCL list size L=4
Decoder config: build/ebn0_sweep.ini (same as the existing static-pm-rank /
trace-learned trials)
New files (this directory), none of which modify any original decoder source:
  - mclass_adaptive_m_dataset.cpp -- dataset generator, reuses
    common/mclass_decoder_wrapper.h (MClassDecoder) unmodified
  - train_adaptive_m.py -- MLP regressor training script (to be added)
CMakeLists.txt was extended (add-only) with a new `mclass_adaptive_m_dataset`
target pointing at this directory; no existing target was changed.

## Idea

Every existing pruning method (static_pm_rank, trace_learned) commits to a
FIXED pruning width k up front, evaluated as its own BLER curve. This
experiment instead trains a model to predict, PER SAMPLE, how wide the
pruning window needs to be -- an adaptive m instead of a fixed one.

Input: the static-pm ranker's own output -- the sorted (ascending) vector of
per-branch checkpoint-t pm_min values across all M branches. This is a
SET-LEVEL feature (dimension M), not a per-branch feature vector: it
deliberately drops which physical branch is which, since the target below is
defined purely by the SHAPE of that sorted score profile (how sharply the
front-runner separates from the pack).

Target (m_required): the smallest m such that keeping the top-m branches,
sorted ascending by checkpoint-t pm_min, is GUARANTEED to contain the branch
that full-PED itself (all M branches, argmin final metric) would have
picked -- and, once included, that branch remains the running best-by-metric
for every larger prefix too. This is well-posed and MONOTONIC: the global
metric minimum is, by construction, also the local minimum of any subset
containing it, so once achieved it cannot be displaced by any branch added
later. m_required is computed by VALUE convergence to the global-minimum
metric (running-min of the prefix reaches full_ped's own metric), not by
matching a single "teacher" branch INDEX -- exact metric ties across
branches are common (independently-converging AED branches that reach the
same codeword report bit-identical metrics), and different tie-break orders
between a natural-index scan and a pm_min-sorted scan can pick different
(but equally valid) tied indices; tracking convergence by value sidesteps
that ambiguity entirely.

## Why this is NOT the same thing as static_bler(k)

This was the single biggest source of confusion while building the dataset
generator, worth recording precisely so it isn't re-litigated.

common/mclass_bler_sim.cpp's static-pm-rank sweep asks, independently for
each FIXED k: "is the running best-by-metric branch among the top-k (by
checkpoint pm_min) correct right now?" This quantity is NOT monotonic in k:
a branch that is the best-so-far AND correct at some small k can later be
DISPLACED by a different, lower-metric but WRONG branch once more branches
enter the prefix (before the true global-metric-minimum branch itself has
entered), making that same k-indexed evaluation flip to incorrect, and it
can flip back to correct again once the actual global optimum finally enters
the prefix. So "smallest k where static_pm_rank(k) happens to be correct"
can be SMALLER than m_required (it can catch an early, transient, lucky hit
that a slightly larger k would have overturned) and is not a well-posed
single threshold per sample.

m_required instead asks a stricter, monotonic question: "smallest m at which
the pruned decode's answer PERMANENTLY locks onto exactly what unrestricted
full decoding would produce." This was confirmed empirically two ways:
  1. A literal re-implementation of mclass_bler_sim.cpp's exact k-loop
     (partial_sort + index ratchet + correct[] check), run inline on the
     same in-memory per-sample decode data as this tool, reproduced
     mclass_bler_sim's real static_bler(k) numbers EXACTLY at every k tested
     (1,2,4,8,16,32,64 -- e.g. 1053/707/423/241/102/41/15 failures out of
     3000 samples at SNR=3.0, checkpoint=32) -- confirming this tool's
     underlying decode data (metrics, checkpoint traces, correctness) is
     byte-for-byte correct, not the source of any earlier mismatch.
  2. That same exact-k-loop replica's implied failure counts (derived from
     m_required alone: fail at k iff m_required>k or full_ped itself is
     wrong) came out systematically HIGHER than the real static_bler(k) at
     every intermediate k -- consistent with the flicker mechanism above,
     since m_required only ever captures the LATER, permanent convergence
     point, not any earlier transient hit.
Confirmed with the user (2026-08-20): m_required (the monotonic, "teacher's
choice survives the pruning" definition) is the intended target, not a
proxy for static_bler(k). The two are expected to diverge and that is not a
bug.

## Practical implication / why this is worth doing

If the regressor predicts m_required exactly, using that many branches
reproduces full_ped's BLER EXACTLY (not just approximately) on every
sample, since by construction pickBest(top-m_required) == full_ped's own
answer. This is a stronger guarantee than any FIXED-k static_pm_rank/
trace_learned curve can offer (those always retain some residual
pruning-induced BLER gap vs. full_ped at any k<M). The research question is
how much AVERAGE compute (mean predicted m, versus M=64) this adaptive
scheme needs to get close to that guarantee, and how that compares to the
average-m needed by fixed-k curves for a similar BLER.

## Dataset generator: mclass_adaptive_m_dataset

    ./mclass_adaptive_m_dataset -ini ebn0_sweep.ini --snr <SNR> \
        --samples <N> --checkpoint <T> \
        --channel-seed <CS> --message-seed <MS> --out <path.csv>

One row per sample:
    sample, M, m_required, full_ped_correct,
    t{T}_pm_min_sorted_1 .. t{T}_pm_min_sorted_M

`full_ped_correct` is kept as a diagnostic column (matches
mclass_bler_sim's full_ped_bler exactly, confirmed at SNR=3.0/checkpoint=32:
15/3000 = 0.005). It is NOT part of m_required's definition (m_required is
defined regardless of whether full_ped itself succeeds -- it's simply the
convergence rank).

## Planned next steps

1. Generate training data at a fixed checkpoint (start with t=32, the
   established sweet-spot-ish middle of the n=128 codeword, mirroring how
   t=8/24 worked out for Golay24) and a training SNR (start with 2.0dB,
   matching the existing trace-learned models' convention), sample size TBD
   based on how much M-dimensional signal is needed -- likely 10-20k
   samples given the input dimension is only M=64.
2. Write train_adaptive_m.py: a small MLP regressor (input dim M=64, a
   couple hidden layers, single scalar output), trained with either plain
   MSE or a loss that penalizes under-prediction more than over-prediction
   (since predicting too LOW risks missing full_ped's answer, while
   predicting too HIGH only costs extra unnecessary compute -- asymmetric
   loss is worth comparing against plain MSE).
3. Evaluate: sweep a safety margin on top of the raw prediction (e.g.
   ceil(pred)+margin, or a trained quantile target) and report, across the
   same SNR grid as the existing sweeps, (a) average m used vs (b) resulting
   BLER (using pickBest over the top predicted-m branches), compared
   against:
     - the fixed-k static_pm_rank/trace_learned curves at the same
       checkpoint (does adaptive beat them at equal average compute?)
     - an "oracle-m" curve using the TRUE m_required per sample (the
       theoretical floor on average-m at zero miss rate -- how close does
       the learned predictor get to it?)
4. Extend to other checkpoints / cumulative-checkpoint features if the
   single-checkpoint version shows a genuine gain, mirroring the
   single-vs-cumulative split already established for Golay24.
