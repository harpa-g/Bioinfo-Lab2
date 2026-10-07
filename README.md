# Oligonucleotide Melting Temperature ($T_m$) Calculator

This repository contains a Python implementation for calculating the melting temperature ($T_m$) of oligonucleotide duplexes based on nucleotide composition.

---

## 1. Calculation for Sequence `ATCGCGTA`

### Sequence Breakdown
* **Sequence ($S$):** `5'-ATCGCGTA-3'`
* **Length ($L$):** 8 base pairs (bp)
* **Base counts:**
  * Adenine ($A$): 2
  * Thymine ($T$): 2
  * Cytosine ($C$): 2
  * Guanine ($G$): 2
* **GC Content:** $\frac{2 + 2}{8} \times 100\% = 50.0\%$

### The Wallace Rule Formula
For short oligonucleotides, the standard heuristic formula is the **Wallace rule** (also known as the Marmur-Schildkraut-Doty short approximation):

$$T_m = 4(G + C) + 2(A + T)$$

Substituting the base counts for $S$:

$$T_m = 4(2 + 2) + 2(2 + 2)$$
$$T_m = 4(4) + 2(4)$$
$$T_m = 16 + 8 = 24\text{ °C}$$

The calculated value is strictly **$24\text{ °C}$**.

---

## 2. Why $T_m$ Cannot Be Higher (~$60\text{ °C}$) for this Sequence

In laboratory PCR protocols, primers are typically designed with an annealing/melting temperature around **$55\text{–}65\text{ °C}$** (often targeting $\sim 60\text{ °C}$). However, sequence `ATCGCGTA` cannot achieve this temperature due to the following biophysical constraints:

### A. Thermodynamic Instability of Short Duplexes (Length Constraint)
The melting temperature represents the temperature at which $50\%$ of the DNA duplex is dissociated into single strands. 

* The stability of a duplex depends directly on the cumulative free energy ($\Delta G = \Delta H - T\Delta S$) contributed by base pairing and base-stacking interactions.
* An 8-mer duplex contains only **7 base-stacking steps** and **20 total hydrogen bonds** ($2 \times 3$ for each GC pair $+ 2 \times 2$ for each AT pair).
* Because total bonding enthalpy ($\Delta H$) scales directly with duplex length, an 8-base sequence does not possess enough cumulative bonding energy to resist thermal denaturation above room temperature ($20\text{–}25\text{ °C}$).

### B. Theoretical Maximum for an 8-mer
Even under the most extreme scenario where an 8-mer is composed of $100\%$ Guanine and Cytosine (e.g., `GCGCGCGC`):

$$T_m^{\text{max}} = 4(8) + 2(0) = 32\text{ °C}$$

Using the Wallace formula, **no 8-nucleotide sequence can mathematically exceed $32\text{ °C}$**, making $\sim 60\text{ °C}$ physically impossible for a duplex of this length under standard conditions.

### C. Required Length for a $60\text{ °C}$ Target
To reach $T_m \approx 60\text{ °C}$ at $50\%$ GC content using the Wallace formula:

$$60 = 4(0.5 \times L) + 2(0.5 \times L)$$
$$60 = 2L + L = 3L$$
$$L = \frac{60}{3} = 20\text{ bp}$$

Standard PCR primers require **18–25 base pairs** to remain stably annealed at typical reaction temperatures ($55\text{–}65\text{ °C}$). An 8-mer would denature immediately under standard PCR annealing conditions.

---

## 3. Summary Table

| Parameter | Sequence `ATCGCGTA` | Typical PCR Primer |
| :--- | :--- | :--- |
| **Length** | 8 bp | 18–25 bp |
| **GC Content** | 50% | 40–60% |
| **Formula Result ($T_m$)** | **$24\text{ °C}$** | **$55\text{–}65\text{ °C}$** |
| **Duplex State at $60\text{ °C}$** | $100\%$ Denatured (Single-stranded) | $\sim 50\%$ Annealed / Duplex |
| **Biological Role** | Micro-oligonucleotide / seed region | Stable PCR primer |

---

## 4. Usage

To run the calculation script:

```bash
python calculate_tm.py
```