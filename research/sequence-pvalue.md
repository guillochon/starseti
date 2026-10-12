# Sequence-match p-value and score bonus

> **Revised.** The exact-substring matcher this note was written for has been retired; `matched.matched_catalog` names every decoded word, in the viewer and in stage 1 alike. The trials and the offset prior below carry over unchanged as the **bits null**; what changed is set out in [What the matched filter changes](#what-the-matched-filter-changes) at the end.

A decoded O–C bit word is named by correlating the folded step means against every window of a packed stream (primes, Fibonacci, counting, powers of two, squares, triangular, each in minimal binary / 4-bit / 8-bit; plus the binary expansions of π and e), over every rotation of the word, at \(12 \le L \le 40\). The word does not have to be exact: a window may disagree with the fold's own two-level word in up to \(\lfloor L/10 \rfloor\) bits and still name it, and the hit records how many. The p-value below is on the **bits**, not the digits; digit count is a display quantity.

## Null

Condition on the decoded word after the MDL bit-period search (do not fold “which \(L\)” into this \(p\)). The \(L\) bits are i.i.d. Bernoulli \(1/2\), independent of the catalogs. Polarity is not searched, so there is no extra \(\times 2\).

> **Revised again.** The fair coin is now conditioned on the word's number of set bits; see [The composition of the word](#the-composition-of-the-word). The trials below are unchanged.

## Trials

The two corpora are different in kind, so they are corrected differently. Half the prior mass goes to each (`SEQ_CORPUS_SPLIT = 2`).

### Numeric packings — flat

A finite corpus, and every alignment in it is a plausible transmission: a beacon sending 11, 13, 17, 19 is no less likely than one sending 2, 3, 5, 7. So the weight is flat over

\[
N(L)=\sum_C \max(0,\,|C|-L+1)
\]

windows, over the 21 packed numeric streams (7 families × 3 encodings). \(N(L)\) is about 4,170. Overlapping windows are positively dependent, so the bound is slightly conservative.

### π and e — weighted by where the match begins

An expansion is **unbounded**, and that breaks a window count in two ways:

1. \(N(L)\) becomes a statement about how many bits we chose to store, not about the hypothesis. The streams used to be 96-bit literals, so every transcendental p-value carried that constant.
2. Storing more punishes every result already obtained. Extending π from 96 to 4,096 bits multiplies a flat count by ~43 and makes every stored hit ~43× less significant — improving the search degrades its own output.

And a window count is the wrong model anyway: **a transmitter sending π starts at the beginning.** Offset 0 is the alignment the hypothesis predicts; a match 74 bits deep is a coincidence with a pseudorandom string, and should be charged as one.

So the offsets carry a decaying prior, \(w_k \propto (k+1)^{-2}\), normalised over an unbounded expansion:

\[
\sum_{k\ge 0}(k+1)^{-2}=\zeta(2)=\frac{\pi^2}{6}
\]

The effective trials for a hit beginning at offset \(k\), over the \(S=3\) expansions searched in parallel, is

\[
T(k)=S\,\zeta(2)\,(k+1)^{2}
\]

The sum converging is the whole point: the expansion can be searched as deep as anyone likes without inflating any p-value, and a hit found 74 bits in pays for having been found 74 bits in.

## p-value

\[
p_{\mathrm{seq}}=\min\!\left(1,\;2\,T\,2^{-L}\right),\qquad
T=\begin{cases}
N(L) & \text{numeric}\\
S\,\zeta(2)\,(k+1)^2 & \text{π / e at offset } k
\end{cases}
\]

Rejecting on any \(p_{\mathrm{seq}}<\alpha\) still controls the family-wise error rate at \(\alpha\) by the union bound (a single-step weighted Bonferroni test; Holm 1979's weighted step-down version would gain little with one hit per star), because the weights sum to at most one. Unmatched, or \(L<6\), or a π/e hit with \(L<12\) → \(p_{\mathrm{seq}}=1\).

| \(L\) | numeric | π/e prefix | π/e from bit 40 | from bit 75 |
| ---: | ---: | ---: | ---: | ---: |
| 12 | 1 | 2.4e-3 | 1 | 1 |
| 16 | 0.13 | 1.5e-4 | 0.24 | 0.85 |
| 20 | 8.0e-3 | 9.4e-6 | 1.5e-2 | 5.3e-2 |
| 24 | 4.9e-4 | 5.9e-7 | 9.4e-4 | 3.3e-3 |
| 32 | 1.8e-6 | 2.3e-9 | 3.7e-6 | 1.3e-5 |

### What changed, and on what evidence

The old formula was a flat union bound over every stream, \(p=\min(1,N(L)2^{-L})\) with \(N(L)\approx 4{,}400\) including the 96-bit expansions. It gave the same p to a prefix match and to one 75 bits deep.

That mattered, because in the store **every** transcendental match is a late-bit match — π at bits 19–80 and e at 51–76 across 87 stars, not one at bit 0, which is the uniform spread over the searched window that chance produces. The two best p-values in the whole store were both of them:

| star | label | \(L\) | old \(p\) | new \(p\) |
|---|---|---:|---:|---:|
| ASASSN 498216382399 | e from bit 75 | 20 | 0.0042 | **0.053** |
| TIC 263054003 | π from bit 40 | 19 | 0.0084 | **0.030** |
| ASASSN 584116059965 | π from bit 2 | 12 | 1.0 | **0.0096** |

Both late hits cleared the 0.01 gate and no longer do; a near-prefix 12-bit hit that the flat bound scored as worthless now does. Since `is_followup` promotes on this gate, that is a change in what reaches stage 2: two chance coincidences out, one near-prefix in.

## The evidence horizon

The expansions are computed rather than stored as literals (`SEQ_STREAM_BITS = 4096`, integer arithmetic, ~5 ms at import; the tests pin the first 96 bits against the literals they replaced). Searching that deep needs a guard: at 4,096 bits **87% of random 12-bit words appear somewhere in π or e**, so a bare `find` would label almost every short word.

The p-value is that guard. A hit is only recorded when \(p_{\mathrm{seq}}<1\), i.e. when \((k+1)^2<2^{L}/(2S\zeta(2))\):

| \(L\) | deepest bit a match can claim | spurious label rate |
| ---: | ---: | ---: |
| 12 | 20 | 0.9% (was 3.8% at 96 bits) |
| 16 | 81 | 0.3% |
| 20 | 325 | 0.0% |
| 24 | 1,303 | 0.0% |
| 32 | 20,860 | 0.0% |

Self-calibrating: the longer the word, the deeper a claim it supports. It is also *stricter* than the old 96-bit literals for short words and far more permissive for long ones, which is the correct trade in both directions. Beyond \(L\approx 26\) the stored depth binds before the prior does; nothing in the data comes close (the longest word in the store is 20 bits).

## Score bonus

Only if the star matched **and** \(p_{\mathrm{seq}} < 0.01\), and only on non-veto classes (not ASTRO / RED NOISE / P-MULTIPLE / SYSTEMATICS / SHAPE / ARTIFACT / ERROR):

\[
e_{\mathrm{seq}}=\min\bigl(1,\; \log_{10}(p_{\mathrm{seq}})/(-6)\bigr),
\qquad \Delta = 0.20\, e_{\mathrm{seq}}
\]

Same log-evidence curve as Ljung-Box / ACF (full scale at \(10^{-6}\)). Added to the non-veto score **before** the sparse-mode penalty. A sparse SURVIVOR still cannot outrank a dense exam.

The same gate is now also a **triage route**: `is_followup` promotes a star whose decoded word clears it, because the autocorrelation statistics are powered by repetition and a word carrying information does not repeat (see [injection-sensitivity.md](injection-sensitivity.md)). Sequence hits still do not mint a classification.

## Stored fields

`stage1_sequence.seq_p` and `seq_bit_len`. The viewer and `score.py` both add \(\Delta\) when `seq_p` is below the gate. Existing rows are backfilled on first `/api/index` (unmatched → \(p=1\); matched → re-decode O–C for \(L\)).

## What this does not cover

- Bits from an O–C mid-threshold are not a fair coin; the p-value is a look-elsewhere correction on the **catalog search**, not a claim that the bits themselves are random.
- Choosing \(L\) via MDL over \(4 \ldots \min(64, n/2)\) is extra look-elsewhere **not** in \(T\). We condition on the decoded word. (`p_decode_shuffle` in stage 1 *does* fold the L search into its null, by searching L on every shuffle draw; the two statistics are complementary.)
- **There is still no across-star correction.** At a 0.01 gate over ~55k stars a few hundred chance hits would be expected if the p-values were uniform; in practice 324 stars match anything at all and one clears the gate, so nothing is on fire, but the gate is not the global statement it looks like. Now that it promotes to triage rather than only adding a score bonus, this is the next thing to fix — a Benjamini–Hochberg pass over the population would say what the false-discovery rate of the promoted set actually is.
- The offset prior's exponent is a **prior**, not a measurement. \((k+1)^{-2}\) is the standard convergent choice and it says a transmitter is far more likely to start at the beginning than 40 bits in. A slower decay (\(1/((k+1)\log^2 k)\)) would be gentler on late hits and is equally valid; the exponent is `SEQ_STREAM_OFFSET_EXPONENT`.
- π is searched both with and without its integer part (a transmitter has no binary point to send), e only with. That asymmetry predates this change and is left alone rather than silently widening the catalog, since another stream raises \(S\) and so every transcendental p-value.

## What the matched filter changes

`matched.matched_catalog` is now the only recognizer. It reads the folded step means, not thresholded bits, so a word whose bits are individually below the noise but collectively right is named, and a long word one or two bits off a catalog window is named as that window with bit errors (the viewer marks the differing bits and says how many). Two nulls apply and the reported \(p\) is the larger of the two p-values:

- **Amplitude.** The unit pattern of a catalog window in the centred folded means, divided by the noise the means carry (total step variance over the count per residue), Bonferroni over every window, rotation and length tried — flat over the numeric packings, by offset over the expansions, times the corpus split. Residuals have to be consistent with that noise (goodness of fit \(\ge 10^{-3}\)) and the window has to fit at least 0.97 as well as the fold's own word. A near miss with \(m\) bits flipped *is* the fold's own word, so its amplitude is that word's. That word is fit to the means, and the best two-level split of \(L\) noise means captures most of their scatter (on pure noise at 33 residues a Gaussian tail charged \(L-1\) ways fell below \(10^{-3}\) nine times in ten), so its amplitude p is the chi-square of the centred means over the noise on \(L-1\) degrees of freedom: the means are not all alike. On noise that sits near 0.4 and clears 0.05 one time in a hundred; TIC 158858696, 33 bits and three off e, went from \(1.1	imes10^{-7}\) to 0.018.
- **Bits.** The formula above, generalised: \(p_{\mathrm{bits}} = \min(1,\; 2\,T\,B(L,m)\,2^{-L})\), with \(B(L,m)=\sum_{i\le m}\binom{L}{i}\) the Hamming ball of the mismatch count. \(T\) counts the distinct numeric windows **at the word's own length** times its \(L\) rotations (a tiled window can only match if its first \(L\) bits do), or \(S\,\zeta(2)\,(k+1)^2\,L\) for an expansion at offset \(k\).

Consequences worth knowing:

- The rotation search is charged. An exact 16-bit numeric hit is not named at all (\(2 \cdot 3{,}000 \cdot 16 \cdot 2^{-16} > 1\)); 17 bits is the shortest numeric word that can carry a name; the 21-bit counting control sits at \(p \approx 0.06\), the 28-bit primes control at \(6.6\times10^{-4}\).
- A 28-bit word one bit off the primes is named at \(p \approx 0.02\), two bits off at \(\approx 0.27\), three bits off not at all. The tolerance is paid for, and a near miss does not reach the 0.01 gate on its own.
- The "evidence horizon" for the expansions is shallower by \(\sqrt{L}\) — for a 20-bit word about 72 bits in rather than 325 — because a rotated word can now be found and is charged for it.
- The word is turned so it reads in catalog order: the fold numbers residues from cycle 0, the matched filter says which rotation lines up with the window, and the per-bit arrays follow.
- `p_catalog_matched` (stage 1) and `seq_p` (viewer, score bonus, triage) are the same statistic computed from two folds of the same steps, and normally agree. Stage 1 has one catalog route: `seq_p` (falling back to `p_catalog_matched` on a payload that predates it) below the 0.01 gate **and** the named window squarer than a sinusoid on the steps; see [The window has to be square](#the-window-has-to-be-square).

## The composition of the word

The bits null above takes the decoded word to be a fair coin: every \(L\)-bit word has probability \(2^{-L}\), so every exact numeric hit of a given length carries the same \(p\), and a 28-bit word of one set bit and 27 clear ones scored 6.6e-4 — the same as the full 28-bit primes control.

That word is not a coincidence with a catalog; it is what one outlier step looks like after the two-level split isolates a single residue. And the powers of two in minimal binary are lone ones between growing runs of zeros: three of their 28-bit windows contain exactly one set bit, and over the 28 rotations searched they cover every position that bit could sit in. **Every 28-bit word with a single set bit is an exact power of two.** ASASSN 429497210103 was named "powers of two [u, 16384]" that way.

The bits from the split are not a fair coin (this note said so under *What this does not cover*), so the null now conditions on what the split produced. Given \(w\) set bits in \(L\), the chance word is uniform over the \(\binom{L}{w}\) arrangements, and the chance that it lies within \(m\) flips of a catalog window with \(w'\) set bits is

\[
\frac{R(L, w', w, m)}{\binom{L}{w}},\qquad
R=\sum_{\substack{a+b\le m\ w'-a+b=w}}\binom{w'}{a}\binom{L-w'}{b}
\]

(\(a\) set bits cleared, \(b\) clear bits set). Summed over the numeric windows at the word's own length, times rotations and the corpus split; for an expansion the same per window, times the offset prior. `matched.reachable` is \(R\), `matched.bits_alignments` the sum, `matched.bits_pvalue` the ratio. Summing \(R\) over every \(w\) recovers the Hamming ball, and summing the alignments over every composition of an exact match recovers the flat count, so this is the same bound with \(\binom{L}{w}\) in place of \(2^L\) and the windows that can be reached counted instead of all of them. On route A (the window itself) the composition is the window's own; on route B (the window with its disagreeing bits flipped, i.e. the fold's own word) it is the fold's word's, per repeat if the window is a tiling.

| word (28 bits) | set bits | fair coin, \(m=0/1/2\) | conditioned, \(m=0/1/2\) |
|---|---:|---|---|
| one spike | 1 | 6.6e-4 / 0.019 / 0.27 | 1 / 1 / 1 |
| two spikes | 2 | 6.6e-4 / 0.019 / 0.27 | 1 / 1 / 1 |
| primes 2–19 | 19 | 6.6e-4 / 0.019 / 0.27 | 1.0e-3 / 0.026 / 0.37 |
| counting 3–9, 4-bit | 13 | 6.6e-4 / 0.019 / 0.27 | 4.7e-4 / 0.015 / 0.21 |
| counting 1–8 (21 bits) | 13 | 0.064 | 0.056 |

Consequences:

- Spike words claim nothing. A word with every bit alike is not a word (\(p=1\)).
- The controls move by less than a factor of two. The 28-bit primes control now sits at **1.0e-3**; the single catalog route gates at 0.01 with the window-squareness guard below (its window reads 0.87), so it reaches triage there as before.
- Balanced words have \(\sim\sqrt{\pi L/2}\) fewer arrangements than \(2^L\), so the expansions' evidence horizon is shallower by that factor again: a 24-bit window of π from 200 bits deep, which the fair coin let through at \(p\approx0.5\), is now past the horizon; one from 60 bits deep is not.
- A near miss at a shallow offset can now outrank an exact hit a few bits deeper. Rotating a window found \(k\) bits into an expansion lines all but \(k\) of its bits up with the prefix, so a shallow exact hit is also a prefix near miss; the prefix costs \((k+1)^2\) less under the offset prior and pays \(R\) for its flipped bits. At \(k=3\), one bit off, the prefix reading wins and is recorded as such (mismatch 1, the bit marked). Under the fair coin the Hamming ball \(\sum_{i\le m}\binom{L}{i}\) always outweighed the prior.

## The window has to be square

Stage 1 used to carry two catalog routes: `p_catalog_matched < 10^{-3}` with `seq_square`, and `seq_p < 0.01` with nothing else. The second was the only promotion route with no two-level test, and it was the one every block-series match took: of the four rows below the gate on 2026-09-19, three (all ASAS-SN, counting words of 26–39 bits) had a null `seq_squareness`, and two of them carried a "coherent secondary period" flag, which is what a smooth wave named by coincidence looks like.

The squareness was null because `_squareness` on the fold's own word needs the word repeated three times — the two-level split was fit to the same folded means it is then judged on, and at fewer repeats it over-explains noise — while the matched filter names a word at two (`y.size >= 2 L k`). Every named word on a short baseline sits in that gap (2.0–2.6 repeats across the matched rows of the store), so requiring `seq_square` on the route would have closed it: 56 of the 1,165 injection recoveries in the sensitivity campaign reached triage on this route alone, all 21-bit counting words at two to three repeats.

The guard is therefore asked of the **catalog window**, not of the fold's own word (`sequence.window_squareness`, `is_square_match`, stored as `window_squareness` in the code and `seq_window_square` in the row). A window is a fixed pattern — the search over windows and rotations is charged in \(p_{\mathrm{seq}}\), not fit here — so two levels against the best sinusoid is a fair contest as soon as the filter has enough of the word to name it. The window is turned back into fold order (window bit \(i\) is fold residue \(i - r \bmod L\)) and scored exactly as the own-word squareness is, against the same 0.25 floor; on an exact match with three or more repeats the two numbers coincide.

| row | repeats | own squareness | window squareness | passes |
| --- | ---: | ---: | ---: | --- |
| ASASSN 575526037701, counting 43–46, 29 bits | 2.5 | — | 0.30 | yes |
| ASASSN 214749425367, counting 27–30, 39 bits | 2.5 | — | 0.00 | no |
| ASASSN 335007454137, counting 6–10, 26 bits | 2.6 | — | 0.17 | no |
| SYNTH cepheid, primes 2–19, 28 bits | 5.3 | 0.87 | 0.87 | yes |
| two-repeat 21-bit counting staircase (test fixture) | 2.0 | — | 0.998 | yes |

The `p_catalog_matched` route is folded into this one: the two statistics are the same matched filter over two folds of the same steps, and the payload's copy is read only when a row predates the stored `seq_p`. `SEQ_CODE_VERSION` is 2; `stage1 --reseq` rewrites the stored codes with the window squareness and `--triage` then re-judges every row from it, so a stored row that predates the change stays out of the catalog route until it is re-decoded rather than being promoted on its old verdict text.

`score.sequence_bonus` still pays on `seq_p < 0.01` alone. It is a ranking bonus on non-veto classes, not a gate, and is left as is.

## Polarity

*Added 2026-09-22.* The filter's sign convention (`_SIGN`: a set bit is a
negative step) is the physics of the trigger, but which of the two symbols a
transmitter kicks on is its choice, not ours: a catalog word sent by a
transmitter that kicks on the clear bits is the same word with every bit
flipped. `matched_catalog` therefore tries each window against the folded
means and against their negative. Everything polarity-free is computed once
(the fold, the noise, the fold's own two-level split and its fit); under the
flipped polarity the means and the free word turn over together, the
window's mismatch to the word becomes L minus what it was, and the free
word's composition becomes L minus its set bits. The best hit carries
`polarity`, `format_label` writes `inverted`, and `sequence.recognize_hit`
complements the window so the viewer marks bits in the word's own
convention; the squareness test is symmetric under the flip.

Every trial count is multiplied by `SEQ_POLARITIES` = 2 alongside the corpus
split (`matched._SPLIT`): the amplitude charge on both routes, the bits
alignments and the stream offset charge. The 28-bit primes control moves
from 1.0e-3 to 2.0e-3 and still clears the 0.01 gate.

Measured before the change on the 12,054 stars at score ≥ 0.5 with a
decodable series (`scratchpad/variants.py`): the upright catalog named 98,
the flipped one 91, and of the 200 unrecognized stars a flipped or reversed
catalog would have named, 5 sat below 0.01 and 1 after charging the search.
Flipping names at the chance rate, which is exactly why it has to be charged
and why it is a freedom rather than a gain. Reversal (the word read
backwards) named at the same chance rate too and has no physical reading --
the fold already searches every rotation and a transmitter sends in time
order -- so it is not tried. Nine candidate extra catalogs (Thue--Morse, the
golden ratio and sqrt 2, a Barker code, an m-sequence, factorials, cubes,
Catalan and Lucas numbers) were tried on the top 3,000 the same day: 13
chance names, the best at p = 0.066, 13 existing names lost to the added
trials. The bit catalogs stay as they are.
