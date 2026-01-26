# SCONE: A General Neural Plug-in Codec for DNA Storage Without Binary Mapping

---

## 📖 Overview

**SCONE** is a trainable and plug-in neural encoder that transforms latent representations (e.g., from VAEs, image compressors) directly into biologically compatible DNA sequences.

Unlike traditional DNA data storage methods which convert binary to quaternary code (bit-to-base mapping), SCONE **bypasses binary entirely**, offering a fully differentiable path from neural representations to DNA with **guaranteed biochemical constraint satisfaction**.

### Key Innovations

- **FSM-Guided Arithmetic Coding**: Deterministic quaternary arithmetic coder with finite state machine (FSM) controlled masking for guaranteed constraint satisfaction
- **Biochemical Constraint Modeling**: Dynamic GC content windowing and homopolymer suppression with fail-safe relaxation
- **Perfect Reversibility**: Lossless encode/decode roundtrip verified over 5000+ sequences

---

## ✨ Key Features

- 📦 **End-to-end differentiable encoding** from neural representations to DNA

- 🧬 **Built-in biochemical constraints** with guaranteed satisfaction:
  - GC content control (windowed, typically 45-55%)
  - Homopolymer suppression (≤3 consecutive identical bases)

- 🔌 **Plug-and-play compatibility** with any upstream model (VAE, COSMOS, CompressAI, etc.)

- 🔁 **Perfect reversibility** via FSM-guided arithmetic coding

- ⚡ **Efficient encoding**: ~1.86 bits/nucleotide with constraint satisfaction

- 🛠️ **Flexible architectures**: Support for MLP and Transformer variants

---

## 🎯 Research Applications

SCONE is designed for research in:

- **Neural image compression & bio-storage**
- **Differentiable biological sequence generation**
- **Next-gen DNA-based archival pipelines**
- **Constraint-aware generative modeling for synthetic biology**

---

## 📂 Project Structure

```
SCONE-DNA/
│
├── config.yaml                    # Training and model configuration
├── requirements.txt               # Python dependencies
├── README.md                      # This file
│
├── src/
│   ├── models/
│   │   └── scone_encoder.py       # Main neural encoder (MLP/Transformer)
│   │
│   └── loss/
│       └── constraint_loss.py     # Biochemical constraint losses
│
├── scripts/
│   └── train_scone.py             # Model training loop
│
├── fsm_constraint.py              # FSM biochemical constraint controller
├── minimal_arithmetic_codec.py    # Standard arithmetic encoder/decoder
├── masked_arithmetic_codec.py     # Masked arithmetic coding with FSM
├── scone_fsm_arith.py             # FSM-guided SCONE arithmetic codec
├── scone_experiments.py           # Camera-ready experiment runner
│
└── experiment_results/            # Saved metrics and CSVs
    ├── metrics.json
    ├── sequence_metrics.csv
    └── summary.csv
```

---

## 🔧 Core Modules

### `FSMConstraint` (fsm_constraint.py)

Deterministic finite state machine for biochemical constraint enforcement.

**Features:**
- **GC Window Control**: Maintains a sliding window (default 20 bases) tracking GC ratio
- **Homopolymer Suppression**: Tracks consecutive base runs, enforces maximum length
- **Fail-safe Relaxation**: Guarantees at least one allowed base by relaxing constraints if mask becomes empty

**Interface:**
```python
fsm = FSMConstraint(gc_window=20, gc_low=0.45, gc_high=0.55, max_homopolymer=3)
mask = fsm.get_mask()      # Returns boolean[4] for allowed bases
fsm.update(base_index)     # Update state after emitting a base
fsm.reset()                # Reset for decoding
```

---

### `MaskedArithmeticCodec` (masked_arithmetic_codec.py)

Arithmetic coding with dynamic probability masking.

**Features:**
- **Masked Renormalization**: Applies boolean mask to probability distribution, renormalizes
- **Frequency Conversion**: Converts probabilities to integer frequencies for arithmetic coding
- **EOS Handling**: Proper end-of-sequence symbol encoding/decoding

---

### `scone_fsm_arith` (scone_fsm_arith.py)

Main FSM-guided SCONE encoder/decoder.

**Interface:**
```python
from scone_fsm_arith import encode_fsm, decode_fsm

# Encode latent symbols to DNA
bitstream, dna_string = encode_fsm(
    latent_symbols,           # List[int] in {0,1,2,3}
    base_probs=[0.25]*4,      # Quaternary probability distribution
    gc_window=20,
    gc_low=0.45, gc_high=0.55,
    max_homopolymer=3
)

# Decode back to latent symbols
decoded_symbols, decoded_dna = decode_fsm(
    bitstream, base_probs, 
    gc_window=20, gc_low=0.45, gc_high=0.55, max_homopolymer=3
)

assert decoded_symbols == latent_symbols  # Perfect reversibility
```

---

### `SCONEEncoder` (src/models/scone_encoder.py)

Neural encoder mapping latent vectors to DNA probability distributions.

**Interface:**
- `forward(latent) → probs`: Generate DNA base probabilities `(batch, seq_len, 4)`
- `sample_sequence(latent) → sequence`: Sample discrete DNA sequence
- `decode_to_string(sequence) → list[str]`: Convert to DNA strings (A/T/G/C)

**Architectures:**
- MLP-based: Fast and lightweight
- Transformer-based: Better sequence dependencies

---

### `Constraint Losses` (src/loss/constraint_loss.py)

Differentiable loss functions for training-time constraint enforcement.

**Available losses:**
1. **`GCContentLoss`**: Controls GC content ratio (typically 40-60%)
2. **`HomopolymerLoss`**: Penalizes consecutive identical bases
3. **`TotalConstraintLoss`**: Weighted combination of multiple constraints

---

## 🚀 Quick Start

### 1. Install Dependencies

```bash
pip install torch numpy matplotlib pyyaml tqdm
# Or install all
pip install -r requirements.txt
```

### 2. Run Experiments

```bash
# Quick test (100 sequences)
python scone_experiments.py --mode quick

# Full experiment (5000 sequences)
python scone_experiments.py --mode full

# Parameter sweep
python scone_experiments.py --mode sweep
```

### 3. Expected Output

```
======================================================================
SCONE EXPERIMENT RESULTS SUMMARY
======================================================================
┌─────────────────────────────────────────────────────────────────┐
│                     SCONE Performance Metrics                    │
├──────────────────────────────┬──────────────────────────────────┤
│ GC Ratio (mean ± std)        │ 0.5002 ± 0.0118                  │
│ Homopolymer (95th pctl)      │                                3 │
│ Bit/nt (mean ± std)          │ 1.8641 ± 0.0688                  │
│ Roundtrip Success Rate       │                          100.00% │
└──────────────────────────────┴──────────────────────────────────┘
```

---

## 📊 Experimental Results

### Large-Scale Evaluation (N=5000, length=100)

| Metric | Value |
|--------|-------|
| GC Ratio | 0.500 ± 0.012 |
| GC Range | [0.46, 0.54] |
| Homopolymer Max | 3 |
| Homopolymer p95 | 3 |
| Bit/nucleotide | 1.864 ± 0.069 |
| Encode Time | 0.60 ms |
| Decode Time | 0.72 ms |
| **Roundtrip Success** | **100%** |

### Parameter Sweep Results

| Config | HP Limit | GC Mean | GC Std | HP p95 | Bit/nt |
|--------|----------|---------|--------|--------|--------|
| Relaxed (0.40-0.60) | 3 | 0.500 | 0.022 | 3 | 2.02 |
| Standard (0.45-0.55) | 3 | 0.500 | 0.012 | 3 | 1.86 |
| Strict (0.48-0.52) | 3 | 0.500 | 0.000 | 3 | 1.28 |

---

## 🧩 Design Philosophy

The system is designed to be **flexible, modular, and plug-and-play** for integration into larger compression-decoding pipelines (e.g., COSMOS, VAE, CompressAI).

**Key principles:**

✅ **Decoupled modules** with clean interfaces

✅ **Configuration-driven** — all hyperparameters configurable

✅ **Guaranteed constraints** — FSM ensures 100% constraint satisfaction

✅ **Perfect reversibility** — lossless encode/decode verified

✅ **Research-ready** — produces publishable metrics and tables

---

## 🔬 Future Extensions

- [ ] Integration with CompressAI learned image compression
- [ ] Structural stability prediction (ViennaRNA, NUPACK)
- [ ] Error-correcting code layers (Reed-Solomon, fountain codes)
- [ ] Sequencing error simulation and robustness evaluation
- [ ] Pre-trained checkpoints for common configurations

---

## 📄 Citation

If you use SCONE in your research, please cite:

```bibtex
@inproceedings{scone2026iscas,
  title={SCONE: FSM-Guided Arithmetic Coding for Constraint-Aware DNA Data Storage},
  author={Author Names},
  booktitle={IEEE International Symposium on Circuits and Systems (ISCAS)},
  year={2026}
}
```

---

## 📧 Contact

For questions or collaboration:
- Open an issue on GitHub
- Email: your.email@example.com

---

## 📜 License

This project is released under the MIT License.

---

> **Note**: This repository contains a fully functional FSM-guided arithmetic coding pipeline with verified 100% roundtrip success over 5000+ sequences.
