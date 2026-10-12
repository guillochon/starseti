# Plan for responding to review 2

*Written 2026-09-25 against `paper/main.tex` at commit 818f6af (the draft the
referee read, plus nothing). Every referee point was checked before a fix was
planned: the physics by recomputation, the literature against the Valency
index, the pipeline claims against `src/` and the `research/` notes, the
counts against the live store (`data/cepheid.sqlite`, which has moved by three
rows since the 23 September snapshot in Table 2: 38,265 flagged, 38,982
vetting rows). No edit has been made. Each item ends with a proposed action;
approve, amend or strike them one at a time. Items marked DECISION need a
call from you before the edit.*

Verdict counts: 3 major, all legitimate; one of them (M2) is confirmed by a
store query to be worse than the referee could see, and one (M3) has a store
answer the referee asked for that turns out to be informative. 9 medium: 8
legitimate, 1 legitimate in conclusion but half-wrong in its argument (m2). 15
minor: all legitimate. Two referee numbers are off (the Ṗ drift in m5 is
three times too large at the Ṗ they quote; the "near half of the chord is
entirely atmosphere" clause of m2 is wrong geometry), neither changes the
conclusion. All five arXiv references the referee cites resolve to the papers
described (Donev & Ivezić 2025, AJ accepted, 228-star table; Netzel et al.
2024, A&A 692, A133; Netzel et al. 2018, MNRAS 480, 1229; Prudil & Skarka
2017, MNRAS; Kovács 2020, ASP Conf. proceedings).

The abstract is at 1919 of 1920 macro-free characters. Every abstract edit
below is paired with a deletion so the total is net-negative.

Two items interact with everything else and should be settled first: the
annular $f_{\rm zone}$ (m2) changes every annular number by a factor of two,
and the gain estimate (M1) may change the reading of the kick ceiling by more
than that. Both are physics, cost no compute, and feed the abstract,
Table 1, Figure 3, Section 6.4 and Conclusions 2 and 3. Do them before M2,
which is the compute-heavy item.

---

## Major

### M1. The gain G is unquantified and the kick ceiling inherits it — LEGITIMATE, load-bearing

**Check.** The referee's scaling is the standard phase-response of a
limit-cycle oscillator to an energy kick: $\Delta\phi/2\pi \approx
Z\,\Delta E/E_{\rm puls}$ with $Z\sim1$ at the best phase, so
$G = Z\,LP/E_{\rm puls}$. $LP$ is $10^{43}$ erg for the reference Cepheid and
$9.4\times10^{39}$ erg for the RRab; both reproduce. The referee's conclusion
("$G=1$ is not obviously conservative") rests on their own rough $E_{\rm
puls}$, which they do not show, and our own two quick estimates disagree with
theirs in the direction of *smaller* $G$:

- Kinetic estimate. With the paper's own $f_{\rm above}\approx10^{-3}$ for
  the RRab and a surface velocity amplitude of $\sim$40 km s$^{-1}$
  ($\sim$65 km s$^{-1}$ peak-to-peak radial velocity times a projection
  factor 1.3, halved), the mass above the driving zone alone carries
  $\tfrac12\,(1.4\times10^{30}\,{\rm g})(4\times10^{6})^2 \approx 10^{43}$ erg,
  a thousand $LP$. That is an upper bound (not all of that mass moves at the
  surface speed) but it is not near unity.
- Energetic bound. In the limit cycle the driving zone's work per cycle
  equals the dissipation, $W = 2\eta E_{\rm puls}$ with $\eta$ the growth
  rate per cycle ($\sim$0.01–0.05 for RR Lyrae fundamental modes), and $W$
  cannot exceed the energy the light curve modulates per cycle,
  $\lesssim 0.3\,LP$. That gives $E_{\rm puls} \lesssim (3$–$15)\,LP$, i.e.
  $G \gtrsim 0.07$–$0.3$ at $Z=1$.

The two bracket $E_{\rm puls}$ between $\sim$3 and $\sim$1000 $LP$ for the
RRab, so $G$ lies somewhere in $10^{-3}$–$0.3$ under the oscillator picture,
and unity is the *optimistic* end. The relaxation-oscillator picture the
paper implicitly uses (the zone fires when its heat content crosses a
threshold, so $\Delta E = L\,\delta t$ advances it by $\delta t$) gives
$G\approx L/L_{\rm dammed}\sim1$–3 instead. Which picture applies is a
question about the phase-response curve the paper already says has never
been computed. The honest statement is a range, and the range straddles
unity from below.

This is the largest physics change in the response, because if $G\ll1$ the
Section 6.4 sentence inverts twice: the TESS RR Lyrae arm would have heard a
self-powered, storage-free transmitter only if $G\gtrsim30$–50 (60–100 after
m2), and the oscillator estimate says $G$ is more likely below unity than
above it. That makes the kick-ceiling statement stronger, not weaker: a
storage-free transmitter at plausible gain is out of reach of any dataset,
existing or planned.

**Cycle-count arithmetic (second half of M1).** At the RRab $G=1$ ceiling,
$\dtp=8\times10^{-5}$ is 3.9 s; against a 20–40 s per-cycle scatter, a
$3\sigma$ residue mean needs $r\approx(3\times30/4)^2\approx500$ repeats of a
28-bit word, $1.4\times10^4$ cycles, 22 yr of continuous per-cycle timing at
0.57 d (the referee's "$\sim$16 yr" uses $10^4$ cycles; same order). Kepler's
4 yr is a factor 4–6 short. Two things the referee did not add and we should:
(i) this assumes the star's clock is white at the 30 s level over decades,
which Section 3.2 says it is not (Blazhko, secular wander), so the real
requirement is worse; (ii) after m2 the ceiling halves and the cycle count
quadruples.

**Proposed actions (M1-a to M1-e):**

- **M1-a.** Section 2.4, after Equation 2: add one paragraph giving both
  estimates above with their inputs (velocity amplitude, $f_{\rm above}$,
  growth rate, modulated luminosity), quoting $G$ as a range and stating that
  $G=1$ is adopted as a reference value, not a conservative one. Cite a
  growth-rate source (Kolláth et al. 2002 or Marconi et al. 2003; verify in
  Valency before adding). DECISION: how deep to go. Recommend the two-estimate
  paragraph (≈120 words) rather than a phase-response calculation, which the
  referee explicitly does not ask for.
- **M1-b.** Label the ceiling "the $G=1$ kick ceiling" at every appearance:
  Section 2.6 (definition), Figure 3 caption and the dash-dotted line label
  in `figures/make_completeness.py`, Section 6.4, Conclusion 3, and the
  viewer's physics box (`src/physics.py`) if it prints the ceiling.
- **M1-c.** Section 6.4: invert the sentence as the referee proposes: "the
  TESS RR Lyrae arm would have heard a self-powered, storage-free transmitter
  only if the star amplifies the beam by $G\gtrsim$ [30–50, or 60–100 after
  m2]", followed by one clause that the estimates of Section 2.4 put $G$ at or
  below unity. Conclusion 3 and the abstract's "30–50" (it is not in the
  abstract; check) follow.
- **M1-d.** Closing paragraph of Section 8: replace "needs more timed cycles
  per star" with the cycle-count arithmetic above (two sentences), so that
  FFI timing and the Gaia DR4 stitch are read as an order-of-magnitude step
  toward a target one to two further orders away unless $G\gg1$.
- **M1-e.** `research/modulation-physics.md`: append a dated section with the
  gain estimates, since the paper says the derivation lives in the notes.

### M2. The store mixes two stage-1 pipelines and the flag rate depends on which — LEGITIMATE, and worse than stated

**Check.** Confirmed by query. The light-curve cache
(`inject.cache_index(survey="tess")`) lists the 1,379 TESS targets that were
fully re-timed with the unclipped timer and RR Lyrae prewhitening
(`research/injection-sensitivity.md`, "Store"); every other TESS row carries
O–C from the earlier timer, with only the routes recomputed from the stored
series. Splitting the dense-mode TESS rows on that membership:

| class, TESS dense | re-timed rows | flagged | not re-timed | flagged | median $N_{\rm cyc}$ (re-timed / not) |
|---|---:|---:|---:|---:|---|
| RR Lyrae | 466 | **54%** | 4182 | 23% | 64 / 77 |
| Cepheid | 435 | 6% | 2267 | 6% | 6 / 16 |
| $\beta$ Cep | 474 | **80%** | 884 | 43% | 510 / 98 |

So the referee's alternative (b), a cycle-count bias in the "any-host" pool,
is ruled out for RR Lyrae: the re-timed stars have *fewer* cycles. It is (a):
Table 2 is out of date for the RR Lyrae and $\beta$ Cep rows, and the 21 per
cent quoted in Sections 4.12 and 6 would be roughly 50 per cent if the whole
sample were timed the way Section 4.3 says it is. The Cepheid arm is
unaffected. The 8,277 TESS RR Lyrae rows in sparse (FFI) mode are a third
population, timed by the block path and never re-timed either. The
"drawn without regard to their stage-1 verdict" phrase is accurate; it is the
comparison with 21 per cent that misleads, because the two numbers come from
different timers.

Also note: the store's TESS RR Lyrae count of 12,925 is 4,648 dense plus
8,277 sparse; the paper never says that most TESS RR Lyrae rows are FFI block
timing, and it should (it changes how the "TESS" column of Table 2 reads, and
it is where the 196 "TESS full-frame" SURVIVORs come from).

**Proposed actions (M2-a to M2-d):**

- **M2-a (preferred; compute).** Re-run stage 1 on the 4,182 un-re-timed
  dense TESS RR Lyrae rows and the 884 $\beta$ Cep rows (the Cepheids can be
  skipped: no change in flag rate, and they contribute nothing to the
  sensitivity). At $\sim$25 s of download plus timing per star, that is
  $\sim$35 h serial, $\sim$9 h at four jobs, and it leaves those light curves
  cached for any later reprocessing. Then re-vet the newly flagged rows with
  stage 2 (the $\sim$1,500 extra TESS RR Lyrae flags and $\sim$300 $\beta$ Cep
  ones; dense vetting is minutes per star). DECISION: run it, or take M2-b.
  Recommend running it, because M2-b leaves the paper's headline counts
  provisional in print and the referee will ask again.
- **M2-b (fallback; text).** Add to Table 2's note which rows were timed
  with which timer, restate the TESS RR Lyrae and $\beta$ Cep flag rates for
  the re-timed subset (54 and 80 per cent) beside the store values, and mark
  the 39 per cent overall and the 11,102 as provisional in Sections 4.12, 6
  and Conclusion 4.
- **M2-c (either way).** Section 4.3: "RR Lyrae light curves from TESS are
  prewhitened before timing" becomes a statement about the current pipeline
  with a clause on which rows carry it; Section 4.11's "92 per cent vetted
  with the revised decoder" gets a parallel sentence about stage 1. After
  M2-a the sentence simplifies to "every TESS RR Lyrae and $\beta$ Cep row".
- **M2-d.** Table 2 or Section 4.4: state the dense/sparse split of the TESS
  rows (4,648 two-minute or 20 s cadence timed per cycle; 8,277 FFI rows
  timed in blocks), so the TESS column is not read as per-cycle throughout.
- **M2-e.** Every downstream count changes after M2-a: Table 2's flagged
  rows, Section 5's store totals, the class tallies of Section 6, the
  11,102, the 21 dense SURVIVORs (Table 4 may change), the 39 per cent, the
  per-survey rates of Section 4.12, the 38,262 in the abstract. Regenerate
  from the store rather than hand-edit; the `util/` check scripts should
  catch stragglers.

### M3. The SURVIVORs have not been tested against Blazhko modulation — LEGITIMATE, and the cheap test is informative

**Check.** The store carries the Gaia DR3 subtype on every target, so the
split the referee asks for is one query. RR Lyrae, current store (flagged =
stage-1 flag; SURVIVOR and ASTRO from stage 2):

| survey | subtype | timed | flagged | SURVIVOR | ASTRO |
|---|---|---:|---:|---:|---:|
| ASAS-SN | RRab | 10,036 | 50% | **9.2%** | 2.9% |
| ASAS-SN | RRc | 4,970 | 40% | **1.0%** | 11.7% |
| ZTF | RRab | 42,192 | 45% | **20.8%** | 2.6% |
| ZTF | RRc | 18,702 | 40% | **5.0%** | 9.8% |
| TESS | RRab | 8,834 | 22% | 2.0% | 2.4% |
| TESS | RRc | 4,013 | 19% | 0.4% | 4.8% |

(RRd rows are too few to matter: 122, 78 and 520.) The answer is
two-sided and better than either of the referee's alternatives. The *flag*
rate is nearly the same for RRab and RRc, so the jitter route is not
measuring Blazhko modulation at the flag stage; it is measuring something
subtype-blind (the linearization bias of Section 4.5, aliasing, Ṗ, or simply
the wrong null). The *survivor* rate is four to nine times higher for RRab
than RRc, in step with the Blazhko incidence ratio (40–50 per cent in RRab
against a few per cent in RRc), and the RRc stars instead land in ASTRO at
three to four times the RRab rate, presumably through the additional-mode
periods the astro-range veto catches. So what survives the sparse vetting is
Blazhko-shaped even though what gets flagged is not. That is exactly the
statement the abstract needs and it costs one table.

The cross-match (referee item 2) is feasible: the Donev & Ivezić table is
228 ZTF field stars with coordinates, and our ZTF sample is the field; the
bulge catalogs (Prudil & Skarka 2017; Skarka et al. 2020) overlap our
$G\le15.5$ sample very little and can be mentioned rather than matched;
Molnár et al. 2022 covers the TESS sample. The match is a cone search the
pipeline already does for TIC and ASAS-SN identifiers.

The coherence test (referee item 3) is a real pipeline addition: the sparse
stage-2 path (`stage2_vet.py`, sparse branch) currently vetoes on secondary
periods of the *block* O–C and on a Keplerian fit to it, but the jitter claim
is vetted only by template-richness, per-camera agreement and season
consistency. A Blazhko star's per-epoch $\delta t_i$ are coherent over the
modulation period; fitting one slow sinusoid (10–1000 d) plus a quadratic to
the $\delta t_i$ and asking whether it explains the jitter variance of
Equation 7 is the sparse analogue of the beat veto. It also absorbs m5 (Ṗ)
in the same fit.

**Proposed actions (M3-a to M3-d):**

- **M3-a (text, no compute).** Add the subtype table above to Section 4.12
  or Section 6 (three survey rows by two subtypes; ten numbers), with the
  two-sided reading: flag rate subtype-blind, survivor rate Blazhko-shaped.
  Cite Jurcsik et al. 2009, Prudil & Skarka 2017, Kovács 2020 for the RRab
  incidence, Netzel et al. 2018 for the RRc incidence, Donev & Ivezić 2025
  for the ZTF lower limit. (All five verified in Valency; add bib entries,
  fields to be checked against ADS per the paper's own acknowledgment.)
- **M3-b (small compute).** Cross-match the SURVIVORs against Donev &
  Ivezić 2025 (ZTF) and Molnár et al. 2022 (TESS); report the matched
  fraction in one sentence. Needs the two tables fetched (VizieR or the
  journal's machine-readable table).
- **M3-c (pipeline; DECISION).** Implement the slow-sinusoid coherence test
  in the sparse stage-2 path: recover $\delta t_i$ per epoch, fit quadratic
  plus one sinusoid over 10–1000 d, veto as ASTRO ("Blazhko-range coherent
  jitter") when the fit explains more than half the excess jitter variance
  and the fitted period is not a sampling alias. Re-vet the 10,885
  ground-based SURVIVORs (ASAS-SN and ZTF light curves are cached; the sparse
  vet is seconds per star, so hours in total). Expected outcome per the
  subtype table: most RRab survivors move to ASTRO. Recommend doing it: it
  is the one change that turns the abstract's headline count into one the
  paper can stand behind, and the referee says so.
- **M3-d (text).** If M3-c is done: abstract and Conclusion 4 report the
  post-test survivor count and say the sparse path now has a Blazhko
  discriminant. If not: abstract's "the vetting cannot attribute to a known
  cause" becomes "the sparse vetting cannot separate from Blazhko modulation"
  (+14 characters; paid for by the m9 abstract edit, −23), and Conclusion 4
  the same.

---

## Medium

### m1. Section 2.2's first clause on the staircase penalty is wrong — LEGITIMATE

**Check.** The running sum of a fair-coin word is a random walk with
excursion $\sim\sqrt{L}\,\delta t/2$ over a word, not "close to its own
ramp"; Section 4.6 says so. The measured factor of two to four comes from
the differencing ($\sqrt2$ in amplitude) and the dropped steps across gaps.

**Action.** Delete "the detrended staircase carries *less* periodic power
than a square wave of the same $\delta t$ would (the running sum of a random
word stays close to its own ramp), and"; keep "each bit now lives in the
difference of two noisy timings". Net shorter.

### m2. Annular $f_{\rm zone}\to1$ is an upper bound by about two — LEGITIMATE in conclusion, half-wrong in argument

**Check.** The referee's first point is right: a ray at $b=r_{\rm zone}$
crosses the zone (thickness $\sim H$, column $\sim\rho_{\rm zone}H$ per
grazing unit) and the overlying atmosphere (column $P/g\approx\rho_{\rm
zone}H$ for a hydrostatic layer) in comparable amounts. The second point is
wrong: the chord's closest approach is at the zone's inner edge, so the zone
occupies the *middle* of the chord and atmosphere both ends; the near half is
atmosphere then zone, not "entirely overlying atmosphere". Redoing the
sequence at $\tau=1$ with equal zone and atmosphere columns: entry atmosphere
absorbs 0.22, zone 0.31, exit atmosphere 0.10, so $f_{\rm zone}=0.31/0.63
\approx0.5$. Aiming the cone $\sim H/2$ deeper puts the closest approach
inside the zone and raises the zone's share toward $\sim$0.7. Adopt
$f_{\rm zone}=0.5$.

**Consequences (DECISION: carry through, or state "of order one half" and
leave the numbers).** Recommend carrying it through, since the referee has
just re-derived all six annular numbers from the paper's inputs and will do
so again. Every annular figure of merit halves:

| quantity | now | at $f_{\rm zone}=0.5$ |
|---|---|---|
| $P_{\rm beam}/L$ annular (Cep 1 EeV / Cep Glashow / RRab / $\beta$ Cep) | 0.09 / $3\times10^{-3}$ / $3.2\times10^{-3}$ / $8.3\times10^{-4}$ | 0.18 / $6\times10^{-3}$ / $6.4\times10^{-3}$ / $1.7\times10^{-3}$ |
| Dyson share at $\eta=0.1$, $\dtp=10^{-3}$ ($\beta$ Cep / RRab / Cep Glashow) | 1 / 3 / 3 per cent | 2 / 6 / 6 per cent; Cepheid at 1 EeV needs more than one sphere |
| Complete-sphere cap ($\beta$ Cep / RRab) | 0.12 / 0.03 | 0.06 / 0.015 |
| $G=1$ kick ceiling ($\beta$ Cep / RRab / Cep) | $3\times10^{-4}$ / $8\times10^{-5}$ / $3\times10^{-6}$ | $1.5\times10^{-4}$ / $4\times10^{-5}$ / $1.5\times10^{-6}$ |
| Section 6.4: erg per bit; $P_{\rm beam}/L$; Dyson share at $\eta=0.1$; ratio to ceiling | $1.5$–$2.5\times10^{38}$; 0.8–1.3$\times10^{-2}$; 8–13%; 30–50 | $3$–$5\times10^{38}$; 1.6–2.6$\times10^{-2}$; 16–26%; 60–100 |
| Annular advantage over central (Cep Glashow / RRab / $\beta$ Cep) | 470 / 260 / 1000 | 230 / 130 / 500 |

The abstract's "15 to 1000 times less power" becomes "a few hundred times
less power" (Glashow taken as the Cepheid optimum, per m9), which is shorter
and drops "depending on the class". Section 6.4's cycle count (M1) quadruples.
Figure 3's dash-dotted lines move (`make_completeness.py`); the viewer's
physics box (`src/physics.py`) changes wherever it carries $f_{\rm zone}$.

**Actions.** Section 2.5: "$f_{\rm zone}\to1$ ... cancelled outright"
becomes "$f_{\rm zone}\approx0.5$ ... reduced from $10^3$ to two"; one clause
that the energy landing above the zone heats the layers that set the surface
luminosity, so a transient light-curve response to each kick is a second
observable. Table 1: add $f_{\rm zone}$ to the note, halve the annular
$P_{\rm beam}/L$ column. Propagate per the table above. Append the
correction to `research/modulation-physics.md`.

### m3. The thermal tail of a kick is not in the injection model — LEGITIMATE

**Check.** `inject.py` offers `staircase` and `square` waveforms only
(`sequence.CODE_MODELS`); every step is instantaneous. The deposited heat
relaxes on the zone's thermal time, of order $P$, so the observable is a
step convolved with a one-to-few-cycle kernel. The fold decoder tolerates
that; the Otsu and squareness statistics are built on clean two-level steps
and would soften.

**Actions.** Section 2.2 or 2.4: one sentence acknowledging the kernel and
calling the clean step an idealization. DECISION on an injection row: a
`--kernel exp:<cycles>` option in `inject.py` that spreads each set bit's
advance over an exponential tail is a small change (the injection is a time
relabelling; the cumulative sum becomes a filtered cumulative sum), and one
TESS RR Lyrae row at a one-cycle tail (40 hosts, the 8-bit word, flagged and
decoded gates) is a few hours of compute. Recommend running it, because the
squareness gate is the discriminant Section 3.3 calls the search's
foundation, and the referee will read "idealization" as "untested".

### m4. Two inconsistencies about vetting triggers and class — LEGITIMATE

**Check.** (a) In `stage1_pulse.py` the catalog route fires on
`seq_p < SEQ_P_GATE` (0.01) *and* `seq_window_square`, i.e. the named
window must beat a free sinusoid on the steps; there is no squareness-free
named-word route. The Table 5 note is therefore right and Section 4.9's "A
named word at $p<0.01$ is also a route to vetting" is incomplete. (b) Stage 2
takes the first hard veto in list order, and the boundary-shift ARTIFACT
test precedes the beat ASTRO test, so TIC 76305643 is ARTIFACT; Section 6.3
item 3 describes both vetoes as if each set the class.

**Actions.** (a) Section 4.9: "A named word at $p<0.01$ whose window beats a
free sinusoid on the steps is also a route to vetting." (b) Section 4.10: one
clause, "the first veto in this order sets the class"; Section 6.3 item 3:
"stage 2 classes the star ARTIFACT on the boundary shift, and would class it
ASTRO on the 2.68-cycle beat had the artifact test passed".

### m5. Ṗ is not removed before the sparse jitter regression — LEGITIMATE

**Check.** `refine_period` in `stage1_pulse.py` solves for a linear period
correction only; the sparse jitter test regresses squared residuals from
that single-period template over the whole baseline, so a secular period
change enters the jitter variance. The referee's number is off: at
$\dot P=10^{-10}$ d d$^{-1}$ over a 5 yr ZTF baseline the quadratic drift is
$\tfrac12(\dot P/P)t^2\approx2.9\times10^{-4}$ d, 25 s, not 80 s; at the more
typical $\dot P\sim3\times10^{-10}$ (0.1 d Myr$^{-1}$) it is 75 s. Either
way it is comparable to the 30–130 s advances the ground-based arm detects at
50 per cent, so the point stands.

**Actions.** Add a $\dot P$ term to `refine_period` (a second regressor,
$g_2 = -\dot F\,x^2/2$, in the same weighted least squares) so the jitter
regression sees a quadratic ephemeris, as the dense path already does; this
is also part of M3-c's fit. Section 4.5: one clause either way. If the
pipeline edit is deferred, list Ṗ beside Blazhko modulation as a known
contributor to the 44–46 per cent flag rate in Section 4.12 and M3-a's
paragraph. Note that Ṗ is subtype-blind, which fits the flag-rate half of
the M3 table.

### m6. Section 7 draws the boundary too generously — LEGITIMATE

**Check.** A compressed code at Table 3 amplitudes adds white timing
variance $(\delta t/2)^2$ per cycle; at 118–197 s that doubles or triples a
quiet TESS RR Lyrae's 20–40 s scatter. A population test of $\sigma_{\rm
O-C}$ against class, period and amplitude was not run. The store has
`oc_rms_seconds`, period and Gaia amplitude on every row, so the test is a
query, not a pipeline change.

**Actions.** Section 7: one sentence, "an efficient message at the same
amplitude would show as excess white jitter against the class's measured
$\sigma_{\rm O-C}$ distribution, a population test this pass did not run; at
the kick ceiling the conclusion stands". Add the population test to the
closing list of next steps. DECISION: optionally run the query now (an hour)
and report the number of $>3\sigma$ outliers in $\sigma_{\rm O-C}$ at fixed
class, period and amplitude for the quiet TESS RR Lyrae; they will mostly
be Blazhko stars and it would give Section 7 a number.

### m7. Gaia DR4 paragraph: one factual error, no source — LEGITIMATE

**Check.** DR4 covers 2014 July 25 to 2020 January 20 (66 months). ZTF began
2018 March; the ASAS-SN $g$ era began 2018. The overlap is about two years,
not "nearly a decade removed".

**Actions.** Replace the clause with "overlaps the first two years of the
ZTF and ASAS-SN $g$ eras and extends four years before them". Cite the ESA
Gaia DR4 release documentation for the 66-month span and the epoch-photometry
projection (add a bib entry of type misc with the URL and access date), or
shrink the paragraph to the two sentences that matter. Recommend the
citation; the paragraph is what the closing recommendation rests on.

### m8. The prewhitening floor is above the modes that contaminate the top ten — LEGITIMATE

**Check.** The floor is 2 per cent of the dominant amplitude; the 0.61 and
0.68 modes are at the per-cent level or below (Netzel et al. 2024 for the
TESS incidence; add to Section 6.1 where the modes are called "the norm").
Both fixes the referee names are pipeline edits: an S/N-based floor in
`prewhiten_modes` (residual-spectrum noise, S/N $\ge4$), and a
common-frequency veto in stage 2 (compare the O–C and per-cycle-amplitude
periodograms the pipeline already computes; Section 6.1 says it "would move
all eight to ASTRO").

**Actions.** DECISION, coupled to M2-a. If M2-a runs, change the floor
first, so the re-timing carries it, and add the common-frequency veto to the
dense stage-2 path before re-vetting; Table 4 is then regenerated and the
text of Section 6.1 describes what the vetting *did*. If M2-a does not run,
implement only the common-frequency veto (stage 2 only, no light curves
needed for the periodogram comparison if the per-cycle amplitude is stored;
check), re-vet the 21 dense SURVIVORs, and rewrite Section 6.1 so Table 4 is
presented as the pre-veto ranking with the veto's outcome stated. Either
way, cite Netzel et al. 2024 ([2404.14373](https://arxiv.org/abs/2404.14373)).

### m9. Two Cepheid cost figures in one paragraph — LEGITIMATE

**Check.** Per bit at each class's optimum, RRab against Cepheid:
$E_{\rm bit}\propto L/(f_{\rm abs}f_{\rm zone})$ gives a ratio of $\sim$1000
only with the Glashow channel for the Cepheid ($f_{\rm abs}=0.15$); with the
1 EeV tunable beam it is $\sim$$3\times10^4$. "A tenth" is the 1 EeV figure.
Conclusion 2 and the abstract's "15 to 1000" mix the same two.

**Actions.** Take the Glashow channel as the Cepheid's optimum everywhere,
since Section 2.5 already says it is the design that applies to the Cepheid:
Section 2.7 "a Cepheid a few $10^{-3}$ with a resonant beam (a tenth with a
tunable one)"; Conclusion 2 the same; abstract "a few hundred times less
power than a central beam" (after m2's factor of two; net −23 characters).

---

## Minor

All checked; one clause each unless noted.

1. Abstract "97 per cent": 10,885/11,102 = 98.0 per cent. Change to 98
   (or restate after M2/M3 change the counts).
2. Section 2.6: $T_{\rm eq}=T_{\rm eff}/\sqrt{86}$ at $43R$ is 740 K for the
   RRab ($T_{\rm eff}=6900$ K) and $\sim$600 K for the Cepheid; write
   "600–750 K". Light-travel anticipation at $43R$: 500 s of 49,000 s is 1.0
   per cent for the RRab, 0.6 per cent for the Cepheid; "well under one per
   cent" becomes "about one per cent or less".
3. Figure 1 caption and Section 2.5: a radial ray through the star crosses
   the shell twice, $2H$, and $2\sqrt{2rH}/2H=\sqrt{2R/H}$; fix "a chord of
   length $H$" in both places.
4. Figure 3: the Cepheid panel has no dash-dotted line because
   $3\times10^{-6}$ (or $1.5\times10^{-6}$ after m2) is off the axis; say so
   in the caption. RR Lyrae panel, 8-bit word at $3.5\times10^{-3}$: decoded
   0.68 against flagged 0.55 (confirmed on the PNG); add "the decode routes
   promote on $p<10^{-3}$ with squareness, so a star can be decoded at 90 per
   cent recall without reaching the flag gate".
5. Table 5 note: the store gives TIC 141590160 at $1.99\times10^{-3}$ and the
   synthetic primes control at $2.01\times10^{-3}$; print both to three
   figures and "would rank ninth" stands.
6. Section 4.1: give the ZTF bright limit (saturation near $r\approx12.5$–13,
   which is why only eight $\beta$ Cep targets have ZTF rows); cite Masci et
   al. 2019 or Bellm et al. 2019 for the limit.
7. Section 4.12: "the empirical false-positive rate" becomes "an upper bound
   on the false-positive rate, exact under the null that no star carries a
   beacon".
8. Section 5: "38,979 vetting results, the excess over the 38,262 flagged
   rows being controls and re-vets" (regenerate both after M2).
9. Section 6.1: add "Table 4 covers ten of the 21 dense-mode SURVIVORs in the
   store, and all ten are mode beating; the dense survivor population is
   essentially one mechanism" (confirmed: 21 TESS dense SURVIVORs). Revisit
   after M2/m8.
10. Section 8 next-steps: add the rule that a parent period longer than a
    third of the baseline defines no harmonic family. DECISION: implement it
    in `stage2_vet.py` (`_harmonic_parents`) now and re-vet AO Tuc, which
    then reaches the beat veto and classes ASTRO; the Section 6.3 text
    changes from "exposes a gap" to "exposed a gap, now closed". Recommend
    implementing; it is a three-line change and it removes the top-ranked
    star's stated failure.
11. `Hippke2014sqrt2` in `references.bib` already says journal "arXiv
    e-prints" with eid 1407.3442; no change. Say so in the response.
12. Acknowledgments: the AI-disclosure paragraph was kept verbatim in
    response 1 as fixed text. The referee has now flagged two factual
    mismatches in it ("equilibrium-solver", "from the solver outputs") and
    "the authors" for a single author. DECISION: (i) replace
    "equilibrium-solver" with "pipeline" and "solver outputs" with "the
    injection database and the store", and "the authors" with "the author"
    throughout the paragraph (minimal, fixes the errors); or (ii) keep
    verbatim and tell the referee it is a fixed disclosure text. Recommend
    (i); a disclosure that names a solver the paper does not have reads as
    boilerplate, which is the opposite of what a disclosure is for.
13. Section 2.4: delete "Transparency is not a fundamental penalty, because
    the cross-section rises with energy, but"; start at "The deep-inelastic
    cross-section is linear in $E$ only to a few TeV".
14. Equation 8: after "$S\,\zeta(2)(k+1)^2L$", add "the inverse of the
    offset weight $w_k$ times the rotations, so that the weighted count is
    what the union bound sees".
15. Section 4.5: typical $K$ from the store (400-star samples): ASAS-SN
    median 19 cycles per block (interquartile 13–27), ZTF 35 (27–49). Add
    "typically 20 cycles for ASAS-SN and 35 for ZTF".

---

## Referee-side corrections to make in the response letter

- m2: the "near half of the chord is entirely overlying atmosphere" clause
  is wrong (the zone sits at the chord's midpoint); the conclusion
  $f_{\rm zone}\approx0.5$ is right for the other reason and we adopt it.
- m5: 80 s is $\sim$25 s at the $\dot P$ quoted; the point stands at a
  typical $\dot P$.
- M2: the cycle-count bias hypothesis (b) is ruled out by the store; it is
  (a), and we say which rows.
- M3: the subtype split does not choose between the referee's two readings;
  it shows both are partly right (flags subtype-blind, survivors
  Blazhko-shaped).
- M1: our own gain estimates land below unity, further from Learned et al.
  than the referee's; we quote a range rather than adopt either.

## Decisions needed (summary)

| item | question | recommendation |
|---|---|---|
| M1-a | depth of the gain estimate | two-estimate paragraph, $G$ quoted as a range |
| m2 | carry $f_{\rm zone}=0.5$ through every annular number | yes |
| M2-a | re-run stage 1 on the 4,182 + 884 un-re-timed TESS RR Lyrae and $\beta$ Cep rows ($\sim$9 h at four jobs) then re-vet | yes |
| M3-c | implement the sparse coherence (Blazhko) veto and re-vet 10,885 ground-based SURVIVORs | yes |
| m3 | add an exponential-tail injection kernel and run one TESS RR Lyrae row | yes |
| m5 | add a $\dot P$ regressor to sparse period refinement | yes (folds into M3-c) |
| m6 | run the population $\sigma_{\rm O-C}$ outlier query | optional |
| m8 | S/N prewhitening floor and common-frequency veto | veto yes; floor only if M2-a runs |
| minor 10 | harmonic-family baseline rule, re-vet AO Tuc | yes |
| minor 12 | edit the fixed AI-disclosure paragraph | yes, minimal |

## Implemented 2026-09-25 (code only; no paper edits)

Every "yes" in the decisions table that is a code change is in the working
tree, formatted, linted and passing the affected tests
(`tests/test_physics.py`, `test_routes.py`, `test_serve_candidates.py`
updated to the halved annular numbers):

- **m2** `src/physics.py`: `F_ZONE_ANNULAR = 0.5` for the annular and
  Glashow geometries, docstrings updated; every coupling, ceiling and viewer
  physics-box figure follows. `paper/figures/make_completeness.py` labels the
  line "$G=1$ kick ceiling" (M1-b).
- **m5** `stage1_pulse.jitter_test(detrend=True)`: a quadratic timing drift
  (period error and Ṗ), linearised through the template slope, is removed
  from the residuals before the jitter regression. Helpers
  `timing_drift_design` and `fit_timing_drift`.
- **M3-c** `stage1_pulse.coherent_jitter_test` and the sparse stage-2 veto
  "ASTRO: per-epoch timing offsets coherent over N d explain X% of the jitter
  variance (Blazhko-range modulation)": one sinusoid over 10–1000 d (capped
  at half the baseline) plus the quadratic, fitted to the per-epoch offsets;
  vetoes when it explains ≥ 50 per cent of the jitter variance, improves
  χ² by ≥ 30 and is not a sampling alias. Guard found by the synthetic
  sparse control: a *repeated* word is itself coherent at the word period
  (the primes control's 28-cycle word was "explained" 75 per cent by a
  15.9 d sinusoid), so the offsets are folded on the fitted period, one bin
  per bit at 1–32 cycles per bit, and the steps are priced as a two-level
  code (ln 2 per bit) against the best wave; a fold that wins by ≥ 5 nats is
  code-like and not vetoed. Checked on synthetic sparse series: Blazhko
  sinusoids at 57 and 120 d and at exactly 28 cycles are vetoed; slow words
  at 8 and 32 cycles per bit win by 75–200 nats; white noise by 0.3. The
  paper's Section 4.10 sparse list gains one item and Section 6's "the sparse
  vetting path has fewer ways" changes.
- **m8** `prewhiten_modes(min_snr=4)` plus a 0.2 per cent absolute floor
  (`PREWHITEN_MIN_SNR`, `PREWHITEN_MIN_FRACTION`), replacing the 2 per cent
  floor for every prewhitened class; modes now record `amplitude_snr`.
  Dense stage 2 gains `common_frequency` (Lomb–Scargle of the per-cycle O–C
  and amplitude, both below FAP 10⁻³ within 0.02 cycles⁻¹) and the veto
  "ASTRO: O-C and per-cycle amplitude both modulated at ν (mode beating)".
- **minor 10** `_harmonic_parents(baseline=...)`: a parent longer than a
  third of the baseline defines no family. AO Tuc's 1.06 d peak now reaches
  the beat veto.
- **m3** `inject --kernel-cycles τ`: each set bit's advance accrues as
  1 − exp(−(n+1)/τ) over the following cycles; rows carry `kernel_cycles`.
- **stage 2 CLI**: `--target-file` (bare names resolve in `data/starlists/`)
  and `--survey`, so the re-vets can be scoped. Target lists written:
  `retime_tess_rrlyrae.txt` (4,182), `retime_tess_bcep.txt` (884),
  `revet_tess_dense_survivors.txt` (111 SURVIVOR + WEAK),
  `revet_ground_survivors.txt` (15,798 ASAS-SN + ZTF SURVIVOR + WEAK).

Not done, by design: M3-b (the Blazhko catalog cross-match needs the two
tables fetched), m6's optional population query, and every text edit.

## Suggested order

1. Physics first (m2, M1, m9, m1, minor 2, 3, 13): text and
   `research/modulation-physics.md` only; fixes the abstract, Table 1,
   Figure 3's lines, Section 6.4, Conclusions 2–3.
2. Pipeline edits before any re-run: $\dot P$ regressor (m5), sparse
   coherence veto (M3-c), common-frequency veto and S/N floor (m8),
   harmonic-family rule (minor 10), injection kernel (m3).
3. Compute: stage-1 re-run on the un-re-timed TESS rows (M2-a); stage-2
   re-vet of everything flagged; the kernel injection row (m3); the Blazhko
   cross-match (M3-b).
4. Regenerate every count, Table 2, Table 4, Table 5's note, Figure 3, and
   the Section 6 tallies from the store; then the remaining text items
   (M2-c/d, M3-a/d, m4, m6, m7, minors).
5. Response letter, with the five referee-side corrections above.

## Results of the runs, and the paper edits made from them (2026-09-26)

- **Re-timing with the S/N prewhitening floor** changed the picture M2 had
  predicted. The 4,182 un-re-timed two-minute TESS RR Lyrae rows flag at 22
  per cent with the noise floor, not the 54 per cent the cached pool showed:
  the 54 per cent was the two per cent floor leaving the 0.61/0.68 modes in.
  The 466 + 474 cached rows were re-timed too (`retime_tess_cached_*.txt`),
  so every two-minute TESS RR Lyrae and β Cep row now carries the noise
  floor; the FFI rows do not. Dense TESS RR Lyrae flag rate 25 per cent.
- **Store now:** 97,725 timed, 77,596 sources, 38,230 flagged; classes among
  flagged rows ARTIFACT 17,313 / SURVIVOR 6,823 / WEAK 5,535 / ASTRO 5,525 /
  RED NOISE 1,981 / SHAPE 509 / SYSTEMATICS 307 / P-MULTIPLE 218 / 19 errors.
  Survivors: ZTF 5,979, ASAS-SN 639, TESS FFI 197, TESS dense 8.
- **Ground re-vet** of the 15,798 first-pass SURVIVOR+WEAK rows: 6,670
  SURVIVOR, 5,520 WEAK, 2,786 ARTIFACT (jitter not reproducing once Ṗ is
  fitted), 820 ASTRO (785 via the coherence test).
- **Common-frequency veto:** 493 ASTRO rows carry it, 32 as sole reason.
  AO Tuc: 0.4866 d mode (ratio 0.685) now prewhitened at 0.34 per cent, O–C
  95 → 42 s, class ARTIFACT with beat, common-frequency and 1.0 d-sinusoid
  reasons behind it. TIC 76305643: 0.1827 d mode (0.615) prewhitened, word
  no longer named (p = 1), ARTIFACT. Top named star is now SX Pup
  (TIC 151932981) at 6.0e-5, ARTIFACT/SHAPE/ASTRO.
- **Injections:** one-cycle thermal tail costs the flag gate ×2.6 (4.0e-3 →
  1.1e-2) and the decode gate ×1.4 (4.2e-3 → 5.8e-3); 6 of 80 controls
  flagged. ZTF run with Ṗ removal + stage 2: raw flag crossings unchanged
  within 10 per cent; of 397 flagged injections stage 2 passed 86 SURVIVOR,
  36 WEAK, 159 ARTIFACT, 113 ASTRO (43 with the coherence test).
- **Paper:** every count, Table 2, new Table 3 row + comments, Table 4
  (six dense + four sparse), the new subtype table, Table 5, Sections 4.3,
  4.9, 4.11, 4.12, 5, 6 (intro, 6.1, 6.2, 6.3 reframed as first pass +
  outcome), Conclusion 4–5 and the closing paragraph updated; abstract
  1842/1920. Figure 5 kept as the first-pass record (caption says so);
  Figure 3 regenerated with the halved ceilings and the G = 1 label. The
  only `\colr{}` left are the three gain numbers of Section 2.4.
- **Not done:** M3-b (Blazhko catalog cross-match), m6's population query,
  minor 12 (fixed disclosure text).

## Section 2 rewritten for deposition at the base of the driving zone (2026-09-26)

The hydrostatic f_above of the first two drafts was wrong by 10^3-10^4; a MESA horizontal-branch profile and two analytic envelopes agree on the corrected values. Section 2, Table 1, the Figure 1 and 3 captions, Sections 2.6-2.7, 4.12, 6.4, 7, Conclusions 2-3, the abstract and physics.py now use the base of the driving zone (thermal time = one period), the Glashow beam for every class, f_zone = 1, efficiencies up to unity and no energy store. Full derivation and numbers: the 2026-09-26 section of modulation-physics.md.
