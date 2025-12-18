import numpy as np
import scipy.fft
import scipy.optimize
from collections import defaultdict
import time
import json
from docopt import docopt

# --- DOCOPT USAGE ---

__doc__ = """LUNA: Latent Untargeted Network Annotation

Usage:
  luna.py --atoms=<atoms_json> --reactions=<reactions_json> --substrates=<substrates_json> --targets=<targets_json> [--ppm=<ppm>] [--max_generations=<max_gen>] [--max_mass=<max_mass>] [--exp_peaks=<peaks_json>] [--bin_res_solver=<bin_res_solver>]

Options:
  --atoms=<atoms_json>          Path to atoms JSON file
  --reactions=<reactions_json>  Path to reactions JSON file (with "common" and "rare" keys)
  --substrates=<substrates_json> Path to substrates JSON file (list of formula dicts, e.g., [{"C":14, "H":9, ...}, ...])
  --targets=<targets_json>      Path to targets JSON file (list of target masses, e.g., [315.0274, ...])
  --ppm=<ppm>                   PPM tolerance for formula search [default: 10.0]
  --max_generations=<max_gen>   Max generations for core network expansion [default: 2]
  --max_mass=<max_mass>         Max mass for candidate generator [default: 550.0]
  --exp_peaks=<peaks_json>      Optional path to experimental peaks JSON (list of [mz, intensity] pairs, e.g., [[315.0274, 0.4], ...])
  --bin_res_solver=<bin_res_solver> Bin resolution for NNLS solver [default: 0.05]
"""

# --- 2. COMPONENT: ISOTOPE MODELER ---

class IsotopeModeler:
    @staticmethod
    def generate_envelope(formula_dict, resolution=0.01, min_abundance=1e-4):
        total_probs = np.array([1.0])
        monoisotopic_mass = sum(ATOM_DATA[el]['mass'] * ct for el, ct in formula_dict.items())
        
        for el, count in formula_dict.items():
            if count == 0: continue
            probs = np.array([x[1] for x in ATOM_DATA[el]['iso']])
            poly_atom = np.poly1d(probs[::-1])
            poly_total = poly_atom ** count
            elem_probs = poly_total.coef
            total_probs = np.convolve(total_probs, elem_probs)
        
        total_probs /= np.sum(total_probs)
        n_peaks = len(total_probs)
        mass_axis = monoisotopic_mass + (np.arange(n_peaks)[::-1] * 1.003355)
        mask = total_probs > min_abundance
        return mass_axis[mask], total_probs[mask]

# --- 3. COMPONENT: CANDIDATE GENERATOR ---

class CandidateGenerator:
    def __init__(self, max_mass=500.0):
        self.MAX_MASS = max_mass
        self.BIN_RES = 0.0001
        self.MAX_PACKED_VAL = 500000
        self.THETA_UNIT = 2 * np.pi / self.MAX_PACKED_VAL
        
        self.LIMITS = {'Na': 1, 'K': 1, 'P': 6, 'S': 6, 'Br': 4, 'Cl': 6, 'F': 30}
        self.spectrum_grid_A = None
        self.spectrum_grid_B = None
        self.config_A = {}
        self.config_B = {}
        self.n_bins = 0
        self._initialize_heteroatom_grids()

    def _build_packing_config(self, atom_order):
        config = {}
        current_mult = 1
        for el in atom_order:
            limit = self.LIMITS[el]
            config[el] = {'mult': current_mult, 'limit': limit}
            current_mult *= (limit + 1)
        return config

    def _generate_grid(self, config):
        accumulator = np.ones(self.n_bins, dtype=np.complex128)
        for el, cfg in config.items():
            mass = ATOM_DATA[el]['mass']
            mult = cfg['mult']
            limit = cfg['limit']
            vec = np.zeros(self.n_bins, dtype=np.complex128)
            k_range = np.arange(limit + 1)
            indices = np.round((k_range * mass) / self.BIN_RES).astype(np.int64)
            valid = indices < self.n_bins
            phases = self.THETA_UNIT * k_range[valid] * mult
            vec[indices[valid]] = np.exp(1j * phases)
            accumulator *= scipy.fft.fft(vec)
        return scipy.fft.ifft(accumulator)

    def _initialize_heteroatom_grids(self):
        print("  [Generator] Pre-calculating Dual Orthogonal FFT Grids...")
        bins_needed = (self.MAX_MASS + 50) / self.BIN_RES
        self.n_bins = int(2**np.ceil(np.log2(bins_needed)))
        atoms = ['Na', 'K', 'P', 'S', 'Br', 'Cl', 'F']
        self.config_A = self._build_packing_config(atoms)
        self.spectrum_grid_A = self._generate_grid(self.config_A)
        self.config_B = self._build_packing_config(reversed(atoms))
        self.spectrum_grid_B = self._generate_grid(self.config_B)
        print(f"  [Generator] Grids Ready. Size: {self.n_bins}")

    def _decode_phase(self, phase_val, config):
        if phase_val < 0: phase_val += 2*np.pi
        packed_int = int(round(phase_val / self.THETA_UNIT))
        counts = {}
        temp = packed_int
        decoders = sorted(config.items(), key=lambda x: x[1]['mult'], reverse=True)
        het_mass = 0.0
        for el, cfg in decoders:
            cnt = temp // cfg['mult']
            temp %= cfg['mult']
            counts[el] = cnt
            het_mass += cnt * ATOM_DATA[el]['mass']
        return counts, het_mass

    def _check_valence(self, counts):
        c = counts.get('C', 0); n = counts.get('N', 0); p = counts.get('P', 0); s = counts.get('S', 0)
        max_capacity = 2 + (2 * c) + (3 * n) + (3 * p) + (4 * s)
        consumers = (counts.get('H', 0) + counts.get('F', 0) + counts.get('Cl', 0) + 
                     counts.get('Br', 0) + counts.get('Na', 0) + counts.get('K', 0))
        return consumers <= max_capacity

    def solve_for_mass(self, target_mass, tolerance_ppm=5.0):
        results = []
        mass_c, mass_n, mass_o, mass_h = (ATOM_DATA[el]['mass'] for el in ['C', 'N', 'O', 'H'])
        max_c = int(target_mass / mass_c)
        
        for c in range(max_c + 1):
            m_c = c * mass_c
            rem_c = target_mass - m_c
            if rem_c < -0.5: break
            max_n = int(rem_c / mass_n)
            for n in range(max_n + 1):
                m_cn = m_c + (n * mass_n)
                rem_cn = target_mass - m_cn
                if rem_cn < -0.5: break
                max_o = int(rem_cn / mass_o)
                for o in range(max_o + 1):
                    m_cno = m_cn + (o * mass_o)
                    rem_cno = target_mass - m_cno
                    if rem_cno < -0.5: break
                    max_h = int((rem_cno + 0.5) / mass_h)
                    for h in range(max_h + 1):
                        m_curr = m_cno + (h * mass_h)
                        gap = target_mass - m_curr
                        if gap < -0.05: break
                        
                        idx = int(round(gap / self.BIN_RES))
                        if 0 <= idx < self.n_bins:
                            for offset in range(-1, 2):
                                check_idx = idx + offset
                                val_A = self.spectrum_grid_A[check_idx]
                                if np.abs(val_A) <= 0.5: continue
                                val_B = self.spectrum_grid_B[check_idx]
                                if np.abs(val_B) <= 0.5: continue
                                
                                counts_A, mass_A = self._decode_phase(np.angle(val_A), self.config_A)
                                counts_B, mass_B = self._decode_phase(np.angle(val_B), self.config_B)
                                
                                if counts_A != counts_B: continue
                                
                                counts = counts_A.copy()
                                counts.update({'C': c, 'H': h, 'N': n, 'O': o})
                                if not self._check_valence(counts): continue

                                total_mass = m_curr + mass_A
                                ppm = ((total_mass - target_mass) / target_mass) * 1e6
                                if abs(ppm) < tolerance_ppm:
                                    counts['ppm'] = ppm
                                    counts['mass'] = total_mass
                                    results.append(counts)
        return results

# --- 4. COMPONENT: REACTION FILTER (Hybrid Phase/Enumeration) ---

class ReactionFilter:
    def __init__(self, substrate_formula, max_generations=2, common_reactions=None, rare_reactions=None):
        self.substrate = substrate_formula
        self.generations = max_generations
        self.BIN_RES = 0.0001
        self.MAX_PACKED_VAL = 500000
        self.THETA_UNIT = 2 * np.pi / self.MAX_PACKED_VAL
        
        # Mapping: Formula_String -> Provenance_List ["Substrate", "+OH", ...]
        self.valid_core_map = {}
        self.core_masses = []
        
        self.common_reactions = common_reactions
        self.rare_reactions = rare_reactions
        
        self.total_grid_A = None
        self.total_grid_B = None
        self.rare_config_A = {}
        self.rare_config_B = {}
        
        self._evolve_common_network()
        self._build_metabolic_grids()

    def _formula_to_str(self, f):
        elements = sorted(f.keys())
        return " ".join([f"{el}{f[el]}" for el in elements if f[el] != 0])
        
    def _calculate_mass(self, f):
        return sum(ATOM_DATA[el]['mass'] * ct for el, ct in f.items())

    def _evolve_common_network(self):
        """Generates the 'Core' metabolic space via BFS."""
        print(f"  [Filter] Evolving Core Network (Common Reactions)...")
        
        sub_str = self._formula_to_str(self.substrate)
        # Queue: (Formula, Path)
        current_gen = [(self.substrate, ["Substrate"])]
        
        self.valid_core_map[sub_str] = ["Substrate"]
        self.core_masses.append(self._calculate_mass(self.substrate))
        
        for g in range(self.generations):
            next_gen = []
            for parent_form, parent_path in current_gen:
                for r_name, delta in self.common_reactions.items():
                    child = parent_form.copy()
                    possible = True
                    for el, change in delta.items():
                        child[el] = child.get(el, 0) + change
                        if child[el] < 0: possible = False; break
                    
                    if possible:
                        f_str = self._formula_to_str(child)
                        if f_str not in self.valid_core_map:
                            new_path = parent_path + [r_name]
                            self.valid_core_map[f_str] = new_path
                            self.core_masses.append(self._calculate_mass(child))
                            next_gen.append((child, new_path))
            current_gen = next_gen
            print(f"    Gen {g+1}: {len(current_gen)} new core metabolites.")
            
        print(f"  [Filter] Core Network Size: {len(self.valid_core_map)} unique formulas.")

    def _build_packing_config(self, keys, limits):
        config = {}
        current_mult = 1
        for k in keys:
            limit = limits[k]['limit']
            config[k] = {'mult': current_mult, 'limit': limit, 'delta': limits[k]['delta']}
            current_mult *= (limit + 1)
        return config

    def _build_metabolic_grids(self):
        """
        Convolves Core Space * Rare Space using FFT.
        Grid A and B used for orthogonal phase encoding of Rare counts.
        """
        print(f"  [Filter] Convolving Core Space with Rare Reactions (Dual Phase)...")
        
        # Calculate grid size
        max_core = max(self.core_masses)
        max_rare = sum(self._calculate_mass(v['delta']) * v['limit'] for v in self.rare_reactions.values())
        max_total = max_core + max_rare + 10.0
        
        n_bins = int(2**np.ceil(np.log2((max_total / self.BIN_RES))))
        
        # 1. Build Core Grid (Sparse Spikes, Real=1.0)
        # We only encode existence of mass here. Formula check happens post-decode.
        core_grid = np.zeros(n_bins, dtype=np.complex128)
        for m in self.core_masses:
            idx = int(round(m / self.BIN_RES))
            core_grid[idx] += 1.0 # Additive for overlapping cores
            
        core_fft = scipy.fft.fft(core_grid)
        
        # 2. Build Rare Grids (Complex Phase = Counts)
        r_keys = list(self.rare_reactions.keys())
        self.rare_config_A = self._build_packing_config(r_keys, self.rare_reactions)
        self.rare_config_B = self._build_packing_config(reversed(r_keys), self.rare_reactions)
        
        def make_rare_fft(config):
            # Start with Identity (Impulse at 0)
            acc_fft = np.ones(n_bins, dtype=np.complex128) # Frequency domain accumulator
            
            for r_name, cfg in config.items():
                delta_mass = self._calculate_mass(cfg['delta'])
                mult = cfg['mult']
                limit = cfg['limit']
                
                vec = np.zeros(n_bins, dtype=np.complex128)
                k_range = np.arange(limit + 1)
                
                indices = np.round((k_range * delta_mass) / self.BIN_RES).astype(np.int64)
                phases = self.THETA_UNIT * k_range * mult
                vec[indices] = np.exp(1j * phases)
                
                acc_fft *= scipy.fft.fft(vec)
            return acc_fft

        rare_fft_A = make_rare_fft(self.rare_config_A)
        rare_fft_B = make_rare_fft(self.rare_config_B)
        
        # 3. Convolution (Multiplication in Freq Domain)
        # Result: Spikes at (Core + Rare) masses, with Phase encoding Rare Counts.
        self.total_grid_A = scipy.fft.ifft(core_fft * rare_fft_A)
        self.total_grid_B = scipy.fft.ifft(core_fft * rare_fft_B)
        self.n_bins = n_bins # Store for indexing

    def _decode_rare_counts(self, phase, config):
        if phase < 0: phase += 2*np.pi
        packed_int = int(round(phase / self.THETA_UNIT))
        
        counts = {}
        temp = packed_int
        decoders = sorted(config.items(), key=lambda x: x[1]['mult'], reverse=True)
        
        for r_name, cfg in decoders:
            cnt = temp // cfg['mult']
            temp %= cfg['mult']
            counts[r_name] = cnt
        return counts

    def check_candidate(self, candidate_dict):
        """
        Validates a candidate using the Pre-Convolved Grid.
        Returns: (True, Provenance_List) or (False, None)
        """
        # 1. Check Mass in Total Grid
        mass = candidate_dict['mass']
        idx = int(round(mass / self.BIN_RES))
        
        if idx >= self.n_bins: return False, None
        
        # Jitter check
        for offset in range(-1, 2):
            check_idx = idx + offset
            val_A = self.total_grid_A[check_idx]
            if np.abs(val_A) < 0.1: continue # Threshold
            
            val_B = self.total_grid_B[check_idx]
            if np.abs(val_B) < 0.1: continue
            
            # 2. Decode Phases
            counts_A = self._decode_rare_counts(np.angle(val_A), self.rare_config_A)
            counts_B = self._decode_rare_counts(np.angle(val_B), self.rare_config_B)
            
            # 3. Collision Check
            if counts_A != counts_B: continue
            
            # 4. Reconstruct Core Formula
            # Candidate - (Rare_Counts * Rare_Delta) = Putative Core
            clean_cand = {k:v for k,v in candidate_dict.items() if k in ATOM_DATA and v != 0}
            putative_core = clean_cand.copy()
            
            rare_path_str = []
            possible = True
            
            for r_name, count in counts_A.items():
                if count > 0:
                    rare_path_str.append(f"{count}x({r_name})")
                    delta = self.rare_reactions[r_name]['delta']
                    for el, change in delta.items():
                        putative_core[el] = putative_core.get(el, 0) - (change * count)
                        if putative_core[el] < 0: possible = False
            
            if not possible: continue
            
            # 5. Validate Core against Map
            core_str = self._formula_to_str(putative_core)
            if core_str in self.valid_core_map:
                # Success!
                core_path = self.valid_core_map[core_str]
                full_path = core_path + rare_path_str
                return True, full_path

        return False, None

# --- 5. COMPONENT: NNLS SOLVER ---

class NNLSSolver:
    def __init__(self, bin_res=0.01):
        self.bin_res = bin_res
    
    def discretize_spectrum(self, masses, intensities, min_mz, max_mz):
        n_bins = int((max_mz - min_mz) / self.bin_res) + 1
        grid = np.zeros(n_bins)
        indices = ((masses - min_mz) / self.bin_res).astype(int)
        valid = (indices >= 0) & (indices < n_bins)
        np.add.at(grid, indices[valid], intensities[valid])
        return grid, n_bins

    def solve(self, experimental_peaks, candidates):
        if not candidates:
            print("  [Solver] No candidates to solve.")
            return []
        exp_mzs = np.array([p[0] for p in experimental_peaks])
        exp_ints = np.array([p[1] for p in experimental_peaks])
        min_mz = min(exp_mzs) - 1.0; max_mz = max(exp_mzs) + 5.0
        b, n_bins = self.discretize_spectrum(exp_mzs, exp_ints, min_mz, max_mz)
        
        A_cols = []; valid_candidates = []
        print(f"  [Solver] Building Design Matrix for {len(candidates)} candidates...")
        for cand in candidates:
            atom_dict = {k:v for k,v in cand.items() if k in ATOM_DATA}
            iso_mzs, iso_ints = IsotopeModeler.generate_envelope(atom_dict)
            col_vec, _ = self.discretize_spectrum(iso_mzs, iso_ints, min_mz, max_mz)
            A_cols.append(col_vec)
            valid_candidates.append(cand)
            
        A = np.vstack(A_cols).T
        print(f"  [Solver] Running NNLS on matrix shape {A.shape}...")
        x, residual = scipy.optimize.nnls(A, b)
        
        results = []
        for i, abundance in enumerate(x):
            if abundance > 0.001:
                res = valid_candidates[i].copy()
                res['abundance'] = abundance
                results.append(res)
        return sorted(results, key=lambda k: k['abundance'], reverse=True)

# --- 6. MAIN EXECUTION ---

def main():
    args = docopt(__doc__)
    
    # Load ATOM_DATA
    with open(args['--atoms'], 'r') as f:
        global ATOM_DATA
        ATOM_DATA = json.load(f)
    
    # Load reactions
    with open(args['--reactions'], 'r') as f:
        reactions_data = json.load(f)
        common_reactions = reactions_data['common']
        rare_reactions = reactions_data['rare']
    
    # Load substrates (list of formula dicts)
    with open(args['--substrates'], 'r') as f:
        substrates = json.load(f)
    
    # Load target masses (list of floats)
    with open(args['--targets'], 'r') as f:
        target_masses = json.load(f)
    
    ppm = float(args['--ppm'])
    max_gen = int(args['--max_generations'])
    max_mass = float(args['--max_mass'])
    bin_res_solver = float(args['--bin_res_solver'])
    
    # Optional experimental peaks
    exp_spectrum = None
    if args['--exp_peaks']:
        with open(args['--exp_peaks'], 'r') as f:
            exp_spectrum = json.load(f)
    
    print("=== Spectral Deconvolution (v38 - Metabolic FFT Convolution) ===")
    
    # Init Generator
    generator = CandidateGenerator(max_mass=max_mass)
    
    # Init ReactionFilters for each substrate
    r_filters = []
    for sub in substrates:
        r_filters.append(ReactionFilter(substrate_formula=sub, max_generations=max_gen,
                                        common_reactions=common_reactions, rare_reactions=rare_reactions))
    
    # Analysis
    print("\n[Step 2] Running Analysis Pipeline...")
    all_candidates = []
    
    for mass in target_masses:
        print(f"\n  >> Analyzing feature at {mass:.4f} Da...")
        
        # A. Chemical Scan
        raw_candidates = generator.solve_for_mass(mass, tolerance_ppm=ppm)
        print(f"     Found {len(raw_candidates)} chemical formulas.")
        
        # B. Biological Filter (FFT Check)
        filtered_candidates = []
        for cand in raw_candidates:
            for r_filter in r_filters:
                is_valid, path = r_filter.check_candidate(cand)
                if is_valid:
                    cand_copy = cand.copy()
                    cand_copy['src_feature_mz'] = mass
                    cand_copy['provenance'] = " -> ".join(path)
                    filtered_candidates.append(cand_copy)
                    break  # Valid for at least one substrate
        
        print(f"     {len(filtered_candidates)} remain after Metabolic Convolution.")
        for fc in filtered_candidates:
             print(f"       Trace: {fc['provenance']}")
        all_candidates.extend(filtered_candidates)
    
    if exp_spectrum:
        # Solver
        print("\n[Step 3] Solving Mixture Abundances...")
        solver = NNLSSolver(bin_res=bin_res_solver)
        solved = solver.solve(exp_spectrum, all_candidates)
        
        print("\n=== FINAL RESULTS ===")
        print(f"{'Formula':<30} | {'Abund':<8} | {'Provenance'}")
        print("-" * 80)
        for res in solved:
            f_str = " ".join([f"{k}{v}" for k,v in res.items() if k in ATOM_DATA and v>0])
            print(f"{f_str:<30} | {res['abundance']:.4f}   | {res['provenance']}")
    else:
        # Output filtered candidates
        print("\n=== FILTERED CANDIDATES ===")
        print(f"{'Formula':<30} | {'MZ':<10} | {'PPM':<8} | {'Provenance'}")
        print("-" * 80)
        for cand in all_candidates:
            f_str = " ".join([f"{k}{v}" for k,v in cand.items() if k in ATOM_DATA and v>0])
            print(f"{f_str:<30} | {cand['src_feature_mz']:<10.4f} | {cand['ppm']:<8.4f} | {cand['provenance']}")

if __name__ == "__main__":
    main()