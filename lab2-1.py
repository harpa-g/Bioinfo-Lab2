def calculate_tm(sequence: str) -> float:
    """Calculate the melting temperature (Tm) of an oligonucleotide

    using the Wallace rule: Tm = 4*(G + C) + 2*(A + T)
    """
    seq = sequence.strip().upper()

    # Validate that only standard DNA nucleotide characters are present
    valid_bases = set("ATCG")
    invalid_chars = set(seq) - valid_bases
    if invalid_chars:
        raise ValueError(
            f"Invalid nucleotide(s) found: {', '.join(invalid_chars)}"
        )

    # Count occurrences of each nucleotide
    count_a = seq.count("A")
    count_t = seq.count("T")
    count_c = seq.count("C")
    count_g = seq.count("G")

    # Wallace rule formula
    tm = 4 * (count_g + count_c) + 2 * (count_a + count_t)

    return tm, {
        "Length": len(seq),
        "A": count_a,
        "T": count_t,
        "C": count_c,
        "G": count_g,
        "GC%": round(((count_g + count_c) / len(seq)) * 100, 2)
        if len(seq) > 0
        else 0,
    }


if __name__ == "__main__":
    # Test sequence
    s = "ATCGCGTA"

    try:
        tm_value, stats = calculate_tm(s)
        print(f"Sequence: {s}")
        print(f"Length:   {stats['Length']} bp")
        print(
            f"Counts:   A={stats['A']}, T={stats['T']}, C={stats['C']}, G={stats['G']}"
        )
        print(f"GC%:      {stats['GC%']}%")
        print(f"Tm:       {tm_value} °C")
    except ValueError as e:
        print(f"Error: {e}")