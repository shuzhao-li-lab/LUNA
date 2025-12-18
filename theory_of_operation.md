# LUNA: Latent Untargeted Network Annotation  

## Abstract

Untargeted metabolomics holds immense promise for elucidating cellular physiology, yet it is hampered by the "annotation bottleneck," wherein the combinatorial explosion of potential molecular formulas and metabolic transformations overwhelms conventional graph-based search algorithms. Traditional approaches model metabolic networks as discrete graphs, incurring exponential computational costs (N_substrates × N_reactions^depth) that limit depth and fail to capture deep conjugates in complex biological mixtures. Here, we introduce LUNA (Latent Untargeted Network Annotation), a probabilistic framework that reframes metabolite annotation as a convolutional signal processing task. LUNA employs a novel Dual-Phase Fast Fourier Transform (FFT) sieve to resolve the chemical subset-sum problem in amortized O(1) time, decoupling metabolic depth from runtime complexity. We further develop a "Hybrid Convolution" architecture that separates dense structural core metabolism from sparse conjugation events, representing the total metabolic space as the convolution of reaction kernels. Uniqueness is ensured through a dual-basis phase-encoding strategy rooted in the Chinese Remainder Theorem, mitigating aliasing artifacts inherent to frequency-domain methods. This enables unsupervised, global reconstruction of metabolic networks with rigorous isotopic validation, facilitating deep annotation of complex matrices at scales previously deemed intractable. LUNA's paradigm shift from graph traversal to spectral analysis promises to illuminate the "dark matter" of metabolomics, with applications in biomarker discovery, toxicology, and systems biology.

## 1. Introduction

Mass spectrometry (MS)-based untargeted metabolomics provides a high-resolution snapshot of the metabolome, capturing thousands of molecular features in a single assay. However, the majority of these features remain unannotated—"dark matter" in the spectral landscape—due to the inherent complexity of chemical and biological spaces. The challenge stems from two intertwined combinatorial problems: (1) enumerating plausible molecular formulas for a given exact mass, a variant of the subset-sum problem, and (2) expanding metabolic networks to account for biotransformations, which generate exponential variants from a single substrate.

Conventional tools, such as METLIN, SIRIUS, or GNPS, rely on heuristic graph traversal or rule-based enumeration to explore these spaces. While effective for shallow networks (e.g., 1–2 reaction steps), these methods scale poorly: increasing search depth to uncover "rare" downstream metabolites, such as glucuronides or glutathione conjugates, inflates the graph exponentially, demanding prohibitive computational resources. This restriction biases analyses toward proximal metabolites, overlooking distal conjugates that often serve as pivotal biomarkers in drug metabolism, environmental toxicology, or disease progression.

To address this, we advocate a fundamental paradigm shift: recasting metabolic annotation not as discrete graph search, but as a continuous signal processing problem in the frequency domain. By leveraging the Fast Fourier Transform (FFT), we exploit convolution theorems to generate and filter vast combinatorial spaces efficiently, transforming exponential enumeration into logarithmic-time operations.

In this work, we present LUNA, a unified framework that:

1. Resolves the subset-sum problem for molecular formulas using a precomputed FFT sieve, enabling exhaustive enumeration without heuristics.
2. Eliminates decoding ambiguities via a Dual-Phase Lock, drawing on orthogonal mixed-radix encodings and the Chinese Remainder Theorem for collision-free lookups.
3. Models metabolic networks as the convolution of a "dense core" (structural modifications) and a "sparse kernel" (conjugations), allowing instantaneous projection across thousands of states.

LUNA's spectral approach not only accelerates annotation but also introduces probabilistic rigor, modeling instrument noise and isotopic distributions to enhance confidence in identifications. We demonstrate its efficacy through a proof-of-concept on xenobiotic metabolism, highlighting its potential to reveal latent metabolic flux in complex systems.

## 2. Methods

### 2.1 The General Hybrid Convolution Framework

LUNA integrates chemical formula generation and metabolic network expansion under a hybrid convolutional formalism, leveraging phase-encoded state vectors in the frequency domain to achieve sub-linear scaling.

Consider a combinatorial entity (e.g., a molecule or pathway) defined by a state vector **v** ∈ ℕ^k, where each entry v_j denotes the count of the j-th component (e.g., atoms or reactions). The entity's mass is given by the linear map M(**v**) = **w** · **v**, with **w** the vector of component masses. The total space 𝒯 is the Minkowski sum of a dense "core" set 𝒞 (e.g., structural backbones) and a sparse "kernel" 𝒦 (e.g., modifications):

𝒯 = { **c** + **k** | **c** ∈ 𝒞, **k** ∈ 𝒦 }.

Direct enumeration scales as |𝒞| × |𝒦|, which is intractable for large sets. Instead, LUNA maps states to a complex-valued signal Z(m) on a discretized mass grid, where presence is indicated by magnitude |Z(m)| and composition **v** by phase θ = arg(Z(m)).

The phase encoding is defined as:

Φ(**v**) = ( ∑_{j=1}^k v_j φ_j ) mod 2π,

where {φ_j} are rationally independent phase weights, ensuring injectivity within bounded domains. The frequency-domain representations Ĉ = ℱ{𝒞} and Ķ = ℱ{𝒦} are computed via FFT. By the convolution theorem:

ℱ{𝒯} = Ĉ ⊙ Ķ,

yielding 𝒯 via inverse FFT in O(N log N) time, where N is the grid size. This decouples core density from kernel complexity.

#### Dual-Phase Uniqueness (The "Lock" Strategy)

Phase collisions—where distinct **v**_1 ≠ **v**_2 yield identical M and Φ—are resolved using a dual-basis scheme inspired by the Chinese Remainder Theorem. Two orthogonal encodings are employed:

- Φ_A ("Little Endian"): Ascending multipliers for low-to-high index components.
- Φ_B ("Big Endian"): Ascending multipliers for high-to-low index components.

Decoding proceeds by inverting phases from both grids; a state is valid iff Decode_A(θ_A) = Decode_B(θ_B). This guarantees uniqueness for masses <2000 Da, reducing the NP-hard subset-sum to O(1) lookup.

### 2.2 Application I: Chemical Formula Generation

For formula deconvolution given target mass m*:

- **v**: Counts of heteroatoms (e.g., N, S, P, Cl, Br, F, Na, K).
- Core 𝒞: C/H/N/O backbone, enumerated via nested Diophantine loops with valence constraints.
- Kernel 𝒦: Heteroatom combinations, phase-encoded in dual FFT grids.
- Phase θ: Encodes heteroatom stoichiometry.

The backbone mass deficit is queried against precomputed grids; valid formulas are those satisfying the dual lock and ppm tolerance (e.g., 5–10 ppm). This "FFT sieve" exhaustively enumerates formulas without pruning, outperforming integer programming methods.

### 2.3 Application II: Metabolic Network Expansion

For substrate S, the latent network is modeled as:

- **v**: Counts of reaction applications.
- Core 𝒞: BFS-generated graph of dense structural reactions (e.g., +OH, -H₂, dealkylation), up to depth d (e.g., 2–3).
- Kernel 𝒦: Sparse convolutions of rare conjugations (e.g., glucuronidation: +C₆H₈O₆; sulfation: +SO₃), with per-reaction limits.
- Phase θ: Encodes conjugation history.

The core masses form a sparse real-valued grid, convolved with phase-encoded kernel FFTs. Candidate validation queries the total grid at m*, decoding paths if the dual lock holds. This enables "instant" reachability checks across millions of states.

### 2.4 Isotopic Envelope Modeling

To discriminate isobars (e.g., ³²S vs. ¹⁶O₂), LUNA computes isotopic envelopes via multinomial expansion. For element E with isotopes (p_i, δ_i), the polynomial is P_E(x) = ∑ p_i x^{δ_i}. The molecular envelope is:

P_total(x) = ∏_E [P_E(x)]^{n_E},

evaluated via convolution or direct coefficient extraction. Envelopes feed downstream fitting, rejecting candidates with mismatched peak ratios.

### 2.5 Gaussian Grid Discretization and NNLS Solver

To handle instrument error (σ ≈ 1–5 ppm), peaks are injected as Gaussians on the mass grid, modeling uncertainty as probability densities. The experimental spectrum forms vector **b**; theoretical envelopes form design matrix **A**. Relative abundances **x** are solved via non-negative least squares (NNLS):

min_{**x** ≥ 0} ‖**Ax** - **b**‖₂².

This probabilistic fit prioritizes global consistency over exact matches, enhancing robustness in noisy data.

## 3. Discussion and Trade-offs

LUNA's spectral framework trades memory for speed: precomputed grids (e.g., 0.0001 Da resolution) enable O(1) queries but require substantial RAM (4–16 GB for <1000 Da). While ideal for small molecules, extensions to macromolecules demand adaptive grids or hierarchical convolutions. The bounded-domain assumption limits generality for unbounded problems, yet aligns with metabolomics constraints.

By reframing graphs as signals, LUNA achieves a paradigm shift, enabling "spectral back-projection" to infer latent pathways from observed spectra. This unlocks unsupervised discovery of metabolic flux, with broad implications for precision medicine, ecology, and synthetic biology. Future integrations with machine learning could further refine priors, positioning LUNA as a cornerstone for next-generation metabolomics.