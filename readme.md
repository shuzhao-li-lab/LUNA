# LUNA: Latent Untargeted Network Annotation

**Reveals Hidden Metabolic Networks via Spectral Back-Projection**

![Python](https://img.shields.io/badge/python-3.8+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-alpha-orange)

## Overview

LUNA is a fast, exhaustive metabolite annotation tool for untargeted metabolomics that combines:

- **Dual-phase FFT sieve** for O(1) molecular formula generation from exact mass
- **Hybrid convolutional metabolic modeling** separating dense core reactions (e.g., hydroxylation, desaturation) from sparse phase-II conjugations (glucuronidation, sulfation, GSH)
- **Collision-free decoding** via orthogonal phase encoding (Chinese Remainder Theorem-inspired)
- **Global mixture deconvolution** using isotopic envelopes and non-negative least squares (NNLS)

By reframing graph-based metabolic network expansion as frequency-domain convolution, LUNA enables deep annotation (multiple generations + rare conjugates) that is computationally intractable with traditional enumeration approaches.

## Key Features

- Exhaustive chemical formula enumeration without heuristics
- Instant metabolic reachability checks across millions of potential biotransformations
- Probabilistic isotopic pattern modeling for discrimination of isobars
- Optional NNLS-based abundance estimation for mixture deconvolution
- Fully configurable via external JSON files (elements, reactions, substrates)

## Installation

```bash
git clone https://github.com/yourusername/luna.git
cd luna

# Recommended: create a virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

pip install numpy scipy docopt
```

Only standard scientific Python packages are required — no heavy dependencies.

## Quick Start

1. Prepare the configuration files (examples provided in `examples/` directory).

2. Run annotation:

```bash
python luna.py \
  --atoms=examples/atoms.json \
  --reactions=examples/reactions.json \
  --substrates=examples/substrates.json \
  --targets=examples/targets.json \
  --ppm=10.0 \
  --max_generations=2 \
  --max_mass=800.0
```

3. (Optional) Include experimental peak list for abundance solving:

```bash
python luna.py \
  ... \
  --exp_peaks=examples/exp_peaks.json \
  --bin_res_solver=0.05
```

### Example Output (Efavirenz demonstration)

```
  >> Analyzing feature at 315.0274 Da...
     Found 87 chemical formulas.
     1 remain after Metabolic Convolution.
       Trace: Substrate

  >> Analyzing feature at 332.0301 Da...
     Found 112 chemical formulas.
     1 remain after Metabolic Convolution.
       Trace: Substrate -> OH_Hydroxylation

  >> Analyzing feature at 491.0595 Da...
     Found 156 chemical formulas.
     1 remain after Metabolic Convolution.
       Trace: Substrate -> Glucuronidation
```

## Input JSON File Formats

All input files are standard JSON. Full working examples are in the `examples/` folder.

### `atoms.json`

Defines elemental monoisotopic masses and isotopic distributions.

**Structure:**
```json
{
  "ElementSymbol": {
    "mass": exact_monoisotopic_mass (float),
    "iso": [
      [mass_offset_1, abundance_1],
      [mass_offset_2, abundance_2],
      ...
    ]
  },
  ...
}
```

- `mass`: the monoisotopic (most abundant isotope) mass used for formula mass calculation
- `iso`: list of [Δmass relative to monoisotopic, natural abundance] pairs
- Abundance values should sum to ≈1.0
- At minimum, include the primary isotope `[0.0, 1.0]` for monoisotopic elements

**Supported elements by default:** C, H, N, O, F, P, S, Cl, Br, Na, K (easily extensible).

### `reactions.json`

Defines metabolic transformations.

**Structure:**
```json
{
  "common": {
    "ReactionName": { "Element": delta_count, ... },
    ...
  },
  "rare": {
    "ReactionName": {
      "delta": { "Element": delta_count, ... },
      "limit": maximum_occurrences_per_molecule (integer)
    },
    ...
  }
}
```

- `"common"`: unlimited applications per generation (used for dense core network via BFS)
  - Delta counts can be positive or negative
- `"rare"`: limited applications (convolved via FFT with phase encoding)
  - Requires both `"delta"` and `"limit"`

### `substrates.json`

List of parent compound molecular formulas from which the metabolic network is expanded.

**Structure:**
```json
[
  { "C": 14, "H": 9, "Cl": 1, "F": 3, "N": 1, "O": 2 },
  { "C": 10, "H": 12, "N": 2, "O": 4 }   // optional second parent
]
```

- Array of dictionaries
- Element symbols as keys, non-negative integer counts as values
- Multiple substrates are supported (network expanded from all simultaneously)

### `targets.json`

List of observed exact masses (m/z values) to annotate.

**Structure:**
```json
[
  315.0274,
  332.0301,
  491.0595
]
```

- Simple array of floats
- Interpreted as neutral monoisotopic masses (adjust manually for adducts if needed)

### `exp_peaks.json` (optional, in development YYMV)

Experimental peak list used for NNLS abundance estimation.

**Structure:**
```json
[
  [mz_1, intensity_1],
  [mz_2, intensity_2],
  ...
]
```

- Array of `[m/z, intensity]` pairs
- Intensities can be arbitrary scale (relative or absolute)
- Only required when performing mixture deconvolution

## Hardware Requirements

LUNA is designed to be memory-efficient:

- **Typical usage (<1000 Da, standard small-molecule elements):** 4–8 GB RAM
- **Extended mass range or many rare conjugations (>1500 Da):** up to 12–16 GB RAM
- **CPU:** Modern multi-core recommended (FFT operations parallelize naturally via NumPy/SciPy)

The precomputed FFT grids scale with mass range and bin resolution (0.0001 Da default). Memory usage is dominated by two complex grids of size ~2²ⁿ (power-of-two for FFT efficiency).

No GPU required — pure CPU NumPy/SciPy.

## Limitations & Known Issues

- Currently optimized for small molecules (<2000 Da)
- Valence checking is basic (Senior's rule approximation)
- Rare reaction limits are global per molecule (not per pathway branch)
- No support yet for charged species or adducts beyond [M+H]+/[M-H]- neutral mass interpretation

## Citation

If you use LUNA in your research, please cite:

> LUNA: Latent Untargeted Network Annotation via Spectral Back-Projection  
> (In preparation — preprint coming soon)

## License

MIT License — feel free to use, modify, and distribute.

## Acknowledgments

Built with inspiration from SIRIUS, GNPS, and the broader metabolomics community. Special thanks to the NumPy/SciPy developers for making high-performance scientific computing accessible.
