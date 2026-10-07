import math
import sys


def calculate_tm_wallace(seq: str) -> dict:
    """Calculates Tm using the Wallace / 2-4 Rule:

    Tm = 4*(G + C) + 2*(A + T)
    Best suited for short oligonucleotides (14-20 bp).
    """
    g = seq.count("G")
    c = seq.count("C")
    a = seq.count("A")
    t = seq.count("T")

    tm = 4 * (g + c) + 2 * (a + t)
    return {"formula": "Wallace Rule", "tm": round(tm, 2), "A": a, "T": t, "C": c, "G": g}


def calculate_tm_salt_adjusted(seq: str, na_conc_molar: float = 0.05) -> dict:
    """Calculates Tm using the salt-adjusted empirical formula:

    Tm = 81.5 + 16.6 * log10([Na+]) + 41 * (%GC) - 600 / length

    Parameters:
    - seq: DNA sequence string (uppercase)
    - na_conc_molar: Monovalent cation concentration [Na+] in mol/L (default 0.05 M = 50 mM)
    """
    length = len(seq)
    if length == 0:
        raise ValueError("Sequence length cannot be 0.")

    g = seq.count("G")
    c = seq.count("C")
    gc_fraction = (g + c) / length

    # log10 of monovalent cation concentration in mol/L
    log_na = math.log10(na_conc_molar)

    tm = 81.5 + 16.6 * log_na + 41 * gc_fraction - (600 / length)
    return {
        "formula": "Salt-Adjusted (Schildkraut-Lifson)",
        "tm": round(tm, 2),
        "na_conc_m": na_conc_molar,
        "gc_percent": round(gc_fraction * 100, 2),
        "length": length,
    }


def analyze_sequence(seq: str, na_conc_molar: float = 0.05):
    clean_seq = seq.strip().upper()
    valid_bases = set("ATCG")
    invalid = set(clean_seq) - valid_bases

    if invalid:
        raise ValueError(f"Sequence contains invalid bases: {', '.join(invalid)}")

    wallace = calculate_tm_wallace(clean_seq)
    salt_adj = calculate_tm_salt_adjusted(clean_seq, na_conc_molar)

    print(f"Sequence:        {clean_seq}")
    print(f"Length:          {len(clean_seq)} bp")
    print(f"Base Counts:     A={wallace['A']}, T={wallace['T']}, C={wallace['C']}, G={wallace['G']}")
    print(f"GC Content:      {salt_adj['gc_percent']}%")
    print("-" * 50)
    print(f"1. Wallace Tm:   {wallace['tm']} °C")
    print(f"2. Salt-Adj Tm:  {salt_adj['tm']} °C (at [Na+] = {na_conc_molar} M)")


if __name__ == "__main__":
    sample_seq = "ACGCGTCGCA"
    analyze_sequence(sample_seq, na_conc_molar=0.05)