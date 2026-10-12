# What the transmitter has to do

*Written 2026-09-17.* The search assumes a modulation model — one-sided phase
advances of order dt/P = 1e-3 — and the viewer converts a measured threshold
into a transmitter power with a coupling constant of 0.1. Neither had been
derived. This note works through the physics behind both, finds one modelling
error in `inject.py`, and puts a number on `TICKLE_COUPLING`.

Throughout, the reference star is a classical Cepheid with M = 5 M☉,
R = 50 R☉, L = 3000 L☉, P = 10 d, against an RRab (0.7, 5, 50, 0.57 d), a
β Cep (10, 7, 15000, 0.20 d) and a δ Sct (1.8, 3, 25, 0.10 d).

## 1. One-sided: the model is right

Learned et al. (2008) trigger the κ-mechanism early and cannot trigger it late:
*"Delivery of a relatively small amount of energy at the right time can thus
trip the cycle early,"* resting on *"near the time of normal triggering there is
increased sensitivity to energy input."*

`inject.inject_phase_code` is `shift = frac·P·bit`, `t' = t − shift`, with
`frac > 0` and `bit ∈ {0,1}`; the built-in control uses the same sign
(`shifts = -SYNTH_INJECT_FRAC * period * code[cyc]` in `stage1_pulse.synthetic`).
Nothing can ever be late. This is correct and needs no change.

Note that the paper's own worked example is a **0.05 cycle** advance on a 2 d
period, and its energy budget (average neutrino power ~1e-6 of stellar output)
implies dt/P ~ 1e-6 at unit gain. `SYNTH_INJECT_FRAC = 1e-3` is commented
"matching the Learned et al. trigger scaling" and is neither of those numbers.
It may be a defensible middle choice; the attribution is wrong.

## 2. Persistent: the model is wrong

Our shift depends only on the *current* cycle's bit, so a set bit fires early
and the next clear bit snaps back to the original ephemeris — a square wave in
O–C. That implies the star remembers an absolute phase it gets restored to.

It does not. A Cepheid is an autonomous nonlinear limit-cycle oscillator, and
phase is the marginal direction along the cycle: amplitude perturbations damp
back onto the limit cycle, phase perturbations have zero Lyapunov exponent and
persist indefinitely. Trigger a pulse early and the clock stays reset; the next
unperturbed pulse arrives one period after the *early* one. The observational
anchor is the Eddington–Plakidis style analysis (1929 for long-period
variables; since applied to Cepheids by Turner, Berdnikov and others): pulsator
O–C behaves as a **random walk built from accumulated period fluctuations**, not
as bounded jitter about a fixed ephemeris.

So the physical O–C is the running sum,

    O−C(n) = −dt · Σ_{k≤n} b_k

a one-sided staircase. Three consequences:

1. ~~**Detection gets easier, so the published null is conservative.**~~ The
   original draft argued that the demeaned running sum of an L-bit word has
   amplitude ~L/4 · dt rather than dt, about 5× more signal for a 20-bit word.
   **Measured, this is wrong for a message.** That excursion belongs to a word
   with long runs (`1111 0000`); a near-fair-coin word keeps its running sum
   close to its own ramp — for the injected `random:8` word the detrended
   staircase has an rms of 0.34 dt against the square wave's 0.48 dt — and the
   bit itself now lives in the difference of two noisy timings, √2 noisier
   than a level. On real TESS RR Lyrae hosts the staircase needs dt/P = 5.5e-3
   to reach triage half the time against 2.5e-3 for the square wave, and
   6.3e-3 to be read against 1.7e-3 ([injection-sensitivity.md](injection-sensitivity.md)).
   The physical waveform is *harder* to find, not easier, unless the message
   is low-entropy; the published square-wave numbers were optimistic by 2–4×.
2. **Decoding breaks.** `sequence.decode_bit_code` is
   `bits = [1 if y < mid else 0]` — it thresholds the *level*. On a staircase
   the level is the integral, not the bits, so `squareness` and
   `recognize_sequence` would never name the word even where Ljung–Box and
   max|ACF| fire hard. The fix is to threshold the **first difference** of O–C.
3. **Detrending eats part of it.** One-sided plus persistent is a monotonic
   ramp whose slope is the duty cycle. `analyse_oc` fits `[1, idx, idx²]` and
   subtracts, so constant duty is absorbed — fine. But *duty changes*, which is
   how a real message varies information density, appear as O–C curvature and
   are absorbed too.

None of this is fixed. It is the main open item this note leaves behind.

## 3. Heating the star is not a constraint

Energy in expands the envelope, lowers the mean density and — since P ∝ ρ̄^(−1/2)
— *lengthens* the period, opposite in sign to the advance being signalled. Using
dP/P = 3ΔE/(GM²/R), a Cepheid at dt/P = 1e-3 and unit gain drifts at +0.25 s/yr;
at the paper's own 0.05-cycle example, +12 s/yr, which sits inside the observed
evolutionary range of ~1–100 s/yr.

There is a clean scale here: at unit gain the time for accumulated heating to
cancel the code's own period offset is **τ_KH/3, independent of amplitude, duty
cycle and period** (offset = d·dt, drift = 3·d·dt/τ_KH). That is ~1700 yr for a
10 d Cepheid, ~20,000 yr for an RRab.

But it is not a limit on the transmitter, because it saturates. After τ_KH the
star simply radiates the injected power (L → L+W) and settles at a permanently
larger radius, total change ΔP/P ≈ 3·d·(dt/P). A civilization signalling for
longer than τ_KH is signalling with a star sitting at a new equilibrium, not
with one that is still drifting. **Recorded and set aside.**

It does matter for us in one way: a constant dP/dt is exactly the shape
`analyse_oc` removes with its quadratic term, so the back-reaction can neither
hide a signal nor fake one.

## 4. The energy cost per bit, and the two penalties

Write the trigger energy as E = L·dt/G, where G is the gain the paper invokes
when it calls the star *"a huge signal amplifier"* and *"as with a spark plug"*
— and then declines to quantify (*"Without detailed stellar modeling, we cannot
make precise estimates"*). Everything below is at G = 1 unless stated.

The beam has to survive two losses. **f_abs** is the fraction absorbed at all;
**f_zone** is the fraction of *that* which lands where it can trigger anything.

### f_abs — transparency, and why it is not fundamental

f_abs = Σ/Λ with Σ = M/πR² the disk-averaged column and Λ = m_p/σ(E) the
interaction column. At 1 TeV, σ ≈ 7e-36 cm² gives Λ = 2.4e11 g/cm², against the
paper's own "1e6 km at solar density" = 1.4e11 g/cm². Their number checks out.

So **f_abs ∝ M/R² ∝ surface gravity**, and scales with beam energy through σ:

| class | log g | f_abs @1 TeV | E for τ=1 | L·P (erg/bit) | bits/yr | rel. cost/bit |
|---|---:|---:|---:|---:|---:|---:|
| classical Cepheid | 1.74 | 1.1e-3 | 912 TeV | 9.9e42 | 37 | 1 |
| **RR Lyrae (RRab)** | **2.89** | **1.5e-2** | **65 TeV** | **9.4e39** | **641** | **7e-5** |
| β Cephei | 3.75 | 1.1e-1 | 9 TeV | 9.9e41 | 1826 | 1e-3 |
| δ Scuti | 3.74 | 1.1e-1 | 9 TeV | 8.3e38 | 3653 | 8e-7 |
| sdB | 5.54 | 1.0 | 0.15 TeV | 2.1e37 | 146k | 2e-9 |
| ZZ Ceti | 7.99 | 1.0 | 0.5 GeV | 2.3e33 | 104k | 3e-13 |

Two readings. There is a hard threshold at **log g ≈ 4.7** where τ = 1 —
everything in the classical instability strip is transparent, compact pulsators
are opaque. And for the transparent ones you simply tune the beam: 912 TeV for a
Cepheid, 65 TeV for an RRab, 9 TeV for a β Cep, all energies nature demonstrably
produces. The paper chose 1 TeV to argue penetration *to the core*, but for
signalling you want τ ≈ 1, not τ ≈ 1e-3. **Transparency is not a real penalty.**

RR Lyrae are nevertheless ~15,000× cheaper per bit than Cepheids at fixed 1 TeV
(14× on coupling × 1053× on L·P), and ~1000× cheaper with each star at its own
optimal energy, since E = L·P·(dt/P) and an RRab is 60× fainter with a 17×
shorter period. Plus 641 bits/yr against 37.

### f_zone — deposition, and why it is the real one

Neutrinos deposit ∝ local density, hence ∝ mass, so a beam through the centre
puts most of its energy in the core and almost none in the He II partial
ionization zone where the κ-mechanism lives. Estimating that zone hydrostatically
— scale height H = kT/(μ m_p g) at T = 4e4 K, column above it Σ = P/g — gives a
mass fraction f_above ≈ 6e-4 for a Cepheid and ~1e-3 for the others.

Normalise the beam to the star's own luminosity:

    P_beam / L  =  duty · (dt/P) / (G · f_abs · f_zone)

At each star's τ=1 energy f_abs ≈ 0.6 **by construction**, so this comes out
essentially **identical for every class** — 1.4 for a Cepheid, 0.83 for the
rest, at G = 1, dt/P = 1e-3, duty 0.5. The class choice buys a large *absolute*
energy saving and buys the *mechanism* nothing.

That is the honest state of the paper's proposal with a central beam: at unit
gain the transmitter is running at roughly the star's own luminosity. That
sounds merely tight, but §5 shows it is worse than tight — a transmitter close
enough to collect the energy it beams cannot gather more than L, so a
central-beam design at unit gain is not underpowered, it is impossible. It
survives only on a gain nobody has computed.

## 5. Shaping the beam: an annulus at the limb

The fix is geometric. Aim a hollow cone — an annulus of impact parameter
b ≈ r_zone — instead of filling the disk. Two things happen, and the second is
the larger.

**The chord lengthens.** A ray tangent to a shell of thickness Δr ≈ H at radius
r traverses 2√(2rH) of it instead of the 2H a radial ray crosses, an enhancement
of √(2R/H):

| class | H/R | √(2R/H) | Σ_disk (g/cm²) | Σ_graze (g/cm²) | E for τ=1 |
|---|---:|---:|---:|---:|---:|
| classical Cepheid | 0.029 | 8.3 | 2.6e8 | 6.5e5 | 2.5 EeV |
| RR Lyrae (RRab) | 0.021 | 9.8 | 3.7e9 | 1.8e7 | 640 PeV |
| **β Cephei** | **0.0020** | **31.5** | 2.7e10 | 4.2e8 | **0.57 PeV** |
| δ Scuti | 0.0048 | 20.4 | 2.6e10 | 2.7e8 | 0.90 PeV |

**The beam never touches the core.** This is the real prize. A ray at impact
parameter b only samples radii ≥ b, so an annulus placed at the driving zone
passes through the zone and the tenuous layers above it and *nothing else*.
f_zone goes from f_above to ≈ 1.

At each geometry's own τ=1 energy the figure of merit f_abs·f_zone therefore goes
from 0.6·f_above to 0.6·1:

> **Net gain = 1/f_above — exactly the deposition penalty, cancelled.**

P_beam/L drops from ~1 to **8.3e-4** at unit gain, for every class. The chord
enhancement is what keeps the price tolerable: the required beam hardness rises
by Σ_disk/Σ_graze = 2/(√(2R/H)·f_above), which is ~400× for a Cepheid rather
than the ~1700× that missing the core would otherwise cost.

**β Cephei is the standout.** Its thin, high-gravity envelope (H/R = 0.002) gives
both the largest chord enhancement and the densest limb, so the annular geometry
works there at **0.57 PeV** — an energy IceCube has actually observed — against
2.5 EeV for a Cepheid. Its Z-bump driving zone also sits deeper (T ~ 2e5 K) and
so plausibly holds more mass, which would help again; that part is a lead, not a
calculation.

### The transmitter is in the system, which makes the geometry easy

The beam has to be aimed from close in, not from across the galaxy: at
P_beam/L ~ 1e-3 the only plausible power source is the star itself, so the
transmitter is a collector array at a standoff of tens of stellar radii. Take
10 AU for a Cepheid and 1 AU for an RRab — both ~43 R★, both with an
equilibrium temperature of ~600–700 K, which a structure can survive.

At that range the angles are enormous and pointing stops being a problem:

| class | standoff | star subtends | cone half-angle | tolerance | anticipation |
|---|---:|---:|---:|---:|---:|
| classical Cepheid | 10 AU (43 R★) | 2.66° | 1.33° | 2.3′ (2.9%) | 4990 s |
| RR Lyrae (RRab) | 1 AU (43 R★) | 2.66° | 1.33° | 1.7′ (2.1%) | 500 s |
| β Cephei | 10 AU (307 R★) | 0.37° | 0.19° | 0.02′ (0.2%) | 4990 s |

Aim a hollow cone of half-angle ~1.3° and hold it to ~3% of itself. That is
arcminutes, not the microarcseconds an interstellar standoff would demand.
The light-travel anticipation is a fixed, known offset of well under 1% of a
cycle, and it is identical for every emitter on a shell of a given radius.

**Tracking.** R varies ~10% over the cycle, moving the target shell by ~3.5 H,
so the cone angle has to chirp in step with the pulsation — which the
transmitter is doing anyway, since the phase is the thing it is modulating.

**Symmetry — and why the swarm solves it for free.** A single limb annulus
illuminates a ring, not a shell. Decomposed about the star's centre that is
l = 0 plus significant l = 2 and higher, so part of the energy would drive
non-radial modes rather than the radial fundamental. But emitters distributed
over 4π each illuminate *their own* limb, and the sum over the array is a
genuinely spherical shell. The collecting area the transmitter needs for power
is the same thing that supplies the angular diversity, so the l ≥ 2
contamination is an artifact of assuming a single beam.

It survives only for a *sparse* array, and then it is **an observational
signature**: excess non-radial mode power phase-locked to the modulation is
what partially-tangential driving would look like from here.

### What powers it — and why the annulus is not optional

The collector has to intercept P_beam/(L·η) of the star's sky, where η is the
efficiency of turning collected starlight into beamed neutrinos:

| geometry | P_beam/L | η = 0.1 | η = 0.01 | η = 1e-3 |
|---|---:|---:|---:|---:|
| **annular**, gain 1 | 8.3e-4 | **0.8%** | **8%** | 83% |
| central, gain 1 | 1.4 | 1400% | 14000% | 1.4e5% |

This is the sharpest result in the note. An annular transmitter at unit gain
needs somewhere between a percent and a tenth of a Dyson sphere — a modest
swarm, and the neutrino production efficiency η is the binding constraint
rather than anything astrophysical. A central-beam transmitter at unit gain
needs to collect **more than the entire output of the star it is driving**,
which it cannot do at any η.

So the annular geometry is not an optimisation. For a transmitter powered by
its own star it is the difference between a modest megastructure and a
thermodynamic impossibility — and it is the reason the unquantified gain G is
load-bearing only for the central-beam case.

## 6. What this changes for the pipeline

**The viewer's coupling constant should be ≥ 1, not 0.1.**
`energyPerBitLsunSec` is `8.64 · L · P_days`, which is
`(dt/P = 1e-3) × (c = 0.1) × 86400` collapsed into one number, while the tooltip
beside it reads "×10 for acoustic coupling" — a factor of 100 apart, or 0.2 in
Sagan K.

The physics says neither. c = 1/(G · f_abs · f_zone) is a *divisor*: you spend
more than L·dt, not less, unless the gain G exceeds 1. At unit gain the
motivated range is **1.7** (annular beam, τ=1) to **2800** (central beam, τ=1).
The committed 0.1 implies G·f_abs·f_zone = 10 — it silently assumes the star
amplifies by an order of magnitude *on top of* perfect coupling. The
energy-per-bit box and every Sagan–Kardashev class derived from it are therefore
optimistic by between 17× and 28,000×.

So the right resolution of the 0.1-versus-×10 discrepancy is that the tooltip is
closer, and both are the wrong side of 1. The constant should carry the geometry
that sets it, and an annular transmitter at unit gain is the defensible
optimistic end.

**Target weighting.** Nothing here argues against Cepheids as the *target*: a
beacon optimises reach, and Cepheids win on both terms of the timing SNR (60×
the luminosity, and a 17× longer period makes dt/P = 1e-3 an 864 s shift rather
than 49 s). The paper's reasoning survives the energy critique. But a
cost-conscious transmitter signalling locally picks RR Lyrae, and — if the
annular geometry is how it is done — β Cep is the class that works at the
softest beam.

That is uncomfortable, because `injection-sensitivity.md` shows the sparse arm
cannot reach dt/P = 1e-3 on an RRab at all: ~470 s of per-epoch timing noise
(~165 s per block) against a 49 s injection, which is why
`SYNTH_INJECT_FRAC_SPARSE` is 10× the dense value. **The class that is ~15,000×
cheapest to modulate is the one we are least sensitive to**, and that belongs in
the sensitivity note as a physics statement rather than only as a scaling choice.

## What this does not cover

- **f_above is model-dependent.** The hydrostatic estimate (H = kT/(μ m_p g) at
  T = 4e4 K, Σ = P/g) is order-of-magnitude. Everything in §4 and §5 scales as
  1/f_above, so a real envelope model would move the numbers, though not the
  structure of the argument. β Cep's Z-bump zone in particular was not computed.
- **The gain G is still unknown.** It is the paper's own gap and this note does
  not close it; it only shows that the annular geometry makes G = 1 sufficient,
  which the central beam does not.
- **The cumulative-O–C model is now implemented** (branch
  `staircase-physics-and-recovery`, 2026-09-17). `inject.inject_phase_code`,
  the built-in controls and `synth_preview` write the running sum;
  `sequence.decode_bit_code` thresholds the first difference by default and
  stage 1's shuffle null is differenced the same way; `physics.py` carries the
  couplings, `P_beam/L` and the self-powered ceiling of §5–6 into the viewer
  and the injection report. The square wave survives as `--waveform square` /
  `model="square"` so the two can be compared on the same hosts. What the
  change did to the measured thresholds is in
  [injection-sensitivity.md](injection-sensitivity.md). §2's third point —
  duty-cycle *changes* absorbed by the quadratic detrend — is still open.

## Correction, 2026-09-24: the cross-section does not stay linear

Review 1 (M1) found that the tau = 1 energies above assume sigma proportional to
E all the way up. The deep-inelastic cross-section is linear only to a few TeV;
above ~10 PeV it goes as sigma ~ 7.8e-36 cm^2 (E/GeV)^0.363 (Gandhi, Quigg, Reno
and Sarcevic 1996, 1998), and the two papers cited for "rises with energy" measure
exactly that sublinear regime. The referee also noticed that the annular energies
in section 5 could not have come from the same law as the central ones: the
ratio Sigma_disk/Sigma_graze = 2/(sqrt(2R/H) f_above) is 401, 205 and 63, so
under any single law 912 TeV -> 366 PeV, 65 TeV -> 13 PeV and 9 TeV -> 570 TeV;
only the beta Cep row matched. Both are errors in this note, now fixed in
`src/physics.py`, which computes everything from the columns:

| path | Sigma (g/cm^2) | tau = 1 (this note) | tau = 1 (Gandhi law) |
|---|---:|---:|---:|
| Cepheid, central | 2.6e8 | 912 TeV | ~100 PeV |
| RRab, central | 3.7e9 | 65 TeV | ~70 TeV |
| beta Cep, central | 2.7e10 | 9 TeV | 9 TeV |
| Cepheid, grazing | 6.5e5 | 2.5 EeV | unreachable (> 1e15 GeV) |
| RRab, grazing | 1.8e7 | 640 PeV | ~2e20 eV |
| beta Cep, grazing | 4.2e8 | 0.57 PeV | ~30 PeV |

Consequences. "f_abs = 0.6 by construction" holds only where tau = 1 is
reachable: every central path (all three energies are inside the spectrum nature
produces, now anchored by KM3NeT's ~220 PeV event of 2025) and the beta Cep
annulus. The RRab and Cepheid annuli are evaluated at 1 EeV, where
sigma = 1.4e-32 cm^2 and f_abs is just the optical depth: 0.16 and 6e-3. So
P_beam/L at dt/P = 1e-3, duty 0.5, unit gain is 8.3e-4 (beta Cep), 3.2e-3
(RRab) and 0.09 (Cepheid), not 8.3e-4 "for every class"; the annulus still
beats the central beam for every class, by 15x rather than 1000x for the
Cepheid. The kick ceiling is per class too: 3.0e-4, 7.7e-5 and 2.8e-6. The
Dyson-sphere shares in section 5's table are the beta Cep row; the RRab needs
3% at eta = 0.1 and the Cepheid most of a sphere.

One physical channel rescues the Cepheid annulus: a monochromatic anti-nu_e beam
at the Glashow resonance, 6.3 PeV, sigma ~ 5e-31 cm^2 per electron, so
Lambda ~ 3.9e6 g/cm^2 in an H/He envelope (0.85 electrons per m_p). On the
grazing Cepheid column tau = 0.17, f_abs ~ 0.15 and P_beam/L ~ 3e-3, a factor 4
from the old number rather than 100. On the RRab (tau ~ 5) and beta Cep
(tau ~ 100) columns it is opaque and would deposit above the zone, so those
classes keep the tunable design. IceCube saw one such event (Nature 591, 220).

Two smaller corrections from the same review. The Kelvin-Helmholtz times in
section 3 were tau_KH/3 mislabelled as tau_KH: GM^2/(RL) is 5.2e3 yr for the
Cepheid and 6.2e4 yr for the RRab, and the cancellation times of 1700 and
2e4 yr are the thirds. And the kick ceiling of the 2026-09-17 addendum in
[injection-sensitivity.md](injection-sensitivity.md) is a per-kick limit, so
no duty cycle enters it; a transmitter able to bank one cycle's intake doubles
it for a fair-coin word.

## Corrections after review 2 (2026-09-25)

**The annular f_zone is one half, not one.** A ray at impact parameter
b = r_zone crosses the driving zone (thickness ~H, column ~rho_zone H along
the grazing chord) in the middle of its path and the hydrostatic atmosphere
above it (column P/g ~ rho_zone H) at both ends. The two columns are
comparable, and at tau ~ 1 the entry atmosphere absorbs first: with equal
zone and atmosphere columns and tau = 1 in total, the entry atmosphere takes
0.22 of the beam, the zone 0.31 and the exit atmosphere 0.10, so
f_zone = 0.31/0.63 ~ 0.5. Aiming the cone ~H/2 deeper raises the zone's
share toward ~0.7. `physics.F_ZONE_ANNULAR = 0.5` now carries this, and every
annular figure of merit halves: P_beam/L at dt/P = 1e-3, duty 0.5, unit gain
is 1.7e-3 (beta Cep), 6.4e-3 (RRab), 0.18 (Cepheid at 1 EeV) and 6e-3
(Cepheid, Glashow); the Dyson shares at eta = 0.1 are 2%, 6% and nearly two
spheres (6% via Glashow); the complete-sphere caps are 0.06 and 0.015; the
G = 1 kick ceilings are 1.5e-4, 3.9e-5 and 1.4e-6. The TESS RR Lyrae
thresholds of 118 and 197 s per bit are then 3e38 and 5e38 erg per bit,
P_beam/L = 1.6-2.6e-2, 16-26% of a sphere at eta = 0.1, a factor 60-100
above the G = 1 ceiling. The referee's second argument for the halving, that
the near half of the chord is entirely atmosphere, is wrong (the zone sits at
the chord's midpoint); the first is enough. The energy deposited above the
zone heats the layers that set the surface luminosity, so a photometric
transient per kick is a second observable.

**The gain G is bracketed, and unity is the optimistic end.** For an
autonomous oscillator an energy kick Delta E at the sensitive phase shifts
the phase by Delta phi / 2 pi ~ Z Delta E / E_puls, Z ~ 1, so
G = Z L P / E_puls: the energy radiated per cycle (1e43 erg Cepheid, 1e40 erg
RRab) over the energy stored in the pulsation. Two rough bounds on E_puls for
the RRab: (i) kinetic, the mass above the zone (f_above ~ 1e-3, 1.4e30 g) at
the surface velocity amplitude (~40 km/s) is ~1e43 erg, an upper bound since
the velocity falls inward; (ii) energetic, in the limit cycle the driving
work per cycle W = 2 eta E_puls with eta ~ 0.01 per cycle (RR Lyrae
hydrodynamic models), and W cannot exceed the modulated luminosity per cycle,
~0.3 L P, so E_puls < (3-15) L P. Together: G ~ 1e-3 to 0.3 at Z = 1. The
relaxation-oscillator picture (the zone fires when its heat content crosses a
threshold, so Delta E = L dt advances it by dt) gives G ~ 1-3 instead. Which
applies is the phase-response curve nobody has computed. The paper now
labels every ceiling "G = 1", and Section 6.4 says the TESS RR Lyrae arm
would have heard a storage-free transmitter only if G > 60-100.

**Reaching the G = 1 ceiling.** dt/P = 4e-5 on the RRab is 2 s against a
20-40 s per-cycle scatter; a 28-bit word read through the fold needs
r ~ (3 x 30 / 2)^2 ~ 2000 repeats for a 3 sigma residue mean, ~6e4 timed
cycles, a century of continuous per-cycle timing of a 0.57 d star with a
white clock. FFI timing and a Gaia DR4 stitch are each one order of magnitude
in cycles toward a target one to two further orders away unless G >> 1.

## The driving zone is shallower than the hydrostatic guess: deposition at its base (2026-09-26)

**The error.** The mass fraction above the He II zone was taken as 6e-4
(Cepheid) and 1e-3 (RRab, beta Cep) from "H = kT/(mu m_p g) at 4e4 K,
Sigma = P/g", with a pressure at the zone that no envelope has: it implies
rho ~ 1e-4 g/cm^3 at 40,000 K in an RRab. Three independent estimates put the
real value three to four orders of magnitude lower:

| star | paper | grey radiative envelope | thermal time = P (Cox) | MESA profile |
|---|---|---|---|---|
| Cepheid, 40 kK | 6e-4 | 2e-5 | 5e-5 | -- |
| RRab, 40 kK | 1e-3 | 3.2e-7 | 3.4e-7 | 1.8e-7 |
| beta Cep, 200 kK | 1e-3 | 1.0e-6 | 0.9e-6 | -- |

The MESA row is Track2_6/profile541.data of Zenodo record 15741530 (Korolik
et al. 2025): a 0.629 Msun, 46.2 Lsun, 6725 K, 5.01 Rsun horizontal-branch
star with 2796 zones, read out of the zip with HTTP range requests. At
40,000 K it has 167 g/cm^2 above, rho = 2.2e-8, H/R = 0.021 (the paper's H/R
was right; its pressure was not); at 55-60 kK, where the thermal time of the
overlying layers is one period, 750 g/cm^2, rho = 6.5e-8, H/R = 0.030,
f_above = 7.6e-7. The grey envelope (Kramers + electron scattering,
T^4 = 3/4 Teff^4 (tau + 2/3)) reproduces it to a factor of two, and so does
the radiative-zero solution P^2 = (2A/8.5) T^8.5.

**The new definition of the target layer.** The base of the driving zone is
the depth at which the thermal time of the overlying layers equals one
period. Above it, heat leaks within a cycle; below it, heat is stored for
many cycles and reaches the zone as a ramp of length ~ that thermal time. It
is therefore the deepest layer at which a kick is a step, and it is the layer
whose overlying energy content is L P, which is what makes G ~ 1 the natural
reference. Adopted values (RRab from MESA, the others from the grey
envelope):

| star | T (K) | Sigma_above (g/cm^2) | f_above | H/R | Sigma_graze = 2 Sigma_above sqrt(2R/H) |
|---|---|---|---|---|---|
| Cepheid | 58,000 | 6.5e3 | 1.0e-4 | 0.042 | 9.0e4 |
| RRab | 55,000 | 7.3e2 | 8e-7 | 0.030 | 1.2e4 |
| beta Cep | 220,000 | 9.3e3 | 1.4e-6 | 0.011 | 2.5e5 |

**Consequences.** No grazing path reaches tau = 1 for any neutrino energy
(1 EeV absorbs 8e-4, 1e-4, 2e-3). The Glashow resonance, Lambda = 3.9e6
g/cm^2, is the design for every class: tau = f_abs = 0.023, 3.0e-3, 0.064.
Because tau << 1 the deposition follows the local density along the chord,
84 per cent within one scale height above the tangent point, i.e. inside the
zone: f_zone = 1 (the review-2 f_zone = 0.5 was for the tau ~ 1 path that no
longer exists). At dt/P = 1e-3, duty 0.5, G = 1: P_beam/L = 0.022, 0.16,
7.8e-3; central beam 8.3, 1040, 600; the annulus wins by 380, 6500 and
77,000. Dyson shares at eta = 1: 2.2%, 16%, 0.8%. The paper now lets eta
approach unity and takes the transmitter to have no energy store, so the
G = 1 kick ceiling is dt/P <= f_abs with a complete sphere: 0.023, 3.0e-3,
0.064. The TESS RRab thresholds (2.4e-3, 4.0e-3) are 0.8 and 1.3 of that
ceiling: E_bit = 7e39 and 1.2e40 erg, P_beam/L = 0.39-0.66, 39-66 per cent of
a sphere at eta = 1. Per bit the RRab is cheapest (3.1e39 erg vs 1.6e40 beta
Cep and 4.3e41 Cepheid); per unit of starlight the beta Cep is.

`physics.py` carries the new f_above and H/R, F_ZONE_ANNULAR = 1,
DEFAULT_GEOMETRY = "glashow", STORAGE_ETA = STORAGE_F_COLL = 1, ETA_GRID =
(1, 0.1, 0.01). Two things were checked and rejected on the way: depositing
below the base buys column only with kernel length (tau_th ~ T^5.25, column
~ T^3.75, so f_abs ~ kernel^0.7), and the prompt mechanical lift of the
envelope by a deep heat pulse is 10-30x weaker than heating the zone at
every depth because the lift is paid against the weight of the overlying
column.
