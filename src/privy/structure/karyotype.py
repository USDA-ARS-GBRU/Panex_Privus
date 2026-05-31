"""k-mer telomere/centromere detection + chromosome structure binning.

Telomere default motifs are the plant repeat ``CCCTAAA`` / its reverse complement
``TTTAGGG`` (Arabidopsis-type), overridable for other systems.  Centromere/satellite
detection finds the densest array of a supplied satellite motif set.  All exact
k-mer matching, pure-Python.

Citations / inspiration: geeViz find_telomeres / find_centromeres; standard plant
telomere repeat (Richards & Ausubel 1988).  See scratch/notes/10_ideas_geeviz.md.
Coordinates are 0-based half-open.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

# Arabidopsis-type plant telomere repeat and its reverse complement.
PLANT_TELOMERE_MOTIFS: tuple[str, ...] = ("CCCTAAA", "TTTAGGG")


@dataclass(frozen=True)
class StructureBin:
    """One chromosome-structure interval."""

    start: int
    end: int
    bin_type: str   # "telomere" | "arm" | "pericentromere" | "centromere"

    @property
    def length(self) -> int:
        return self.end - self.start


def _find_motif_positions(sequence: str, motif: str) -> list[tuple[int, int]]:
    spans: list[tuple[int, int]] = []
    if not motif:
        return spans
    start = 0
    while True:
        i = sequence.find(motif, start)
        if i < 0:
            break
        spans.append((i, i + len(motif)))
        start = i + 1
    return spans


def kmer_blocks(
    sequence: str,
    motifs: Sequence[str],
    *,
    max_gap: int = 100,
) -> list[tuple[int, int]]:
    """Merge exact occurrences of any *motif* into blocks (gap ≤ *max_gap*).

    Returns non-overlapping ``(start, end)`` blocks sorted by start.
    """
    spans: list[tuple[int, int]] = []
    for motif in motifs:
        spans.extend(_find_motif_positions(sequence, motif))
    if not spans:
        return []
    spans.sort()
    blocks: list[tuple[int, int]] = []
    cur_start, cur_end = spans[0]
    for s, e in spans[1:]:
        if s - cur_end <= max_gap:
            cur_end = max(cur_end, e)
        else:
            blocks.append((cur_start, cur_end))
            cur_start, cur_end = s, e
    blocks.append((cur_start, cur_end))
    return blocks


@dataclass(frozen=True)
class TelomereResult:
    """Telomere repeat blocks and whether each chromosome end is capped."""

    blocks: list[tuple[int, int]]
    five_prime_capped: bool
    three_prime_capped: bool


def find_telomeres(
    sequence: str,
    motifs: Sequence[str] = PLANT_TELOMERE_MOTIFS,
    *,
    end_window: int = 10_000,
    min_array_len: int = 100,
    max_gap: int = 100,
) -> TelomereResult:
    """Detect telomeric repeat arrays and whether either sequence end is capped.

    An end is "capped" when a repeat block of at least *min_array_len* bp lies
    within *end_window* bp of that end.
    """
    n = len(sequence)
    blocks = [b for b in kmer_blocks(sequence, motifs, max_gap=max_gap)
              if b[1] - b[0] >= min_array_len]
    five = any(start < end_window for start, _end in blocks)
    three = any(end > n - end_window for _start, end in blocks)
    return TelomereResult(blocks=blocks, five_prime_capped=five, three_prime_capped=three)


def find_centromere(
    sequence: str,
    satellite_motifs: Sequence[str],
    *,
    max_gap: int = 1000,
    min_array_len: int = 1000,
) -> tuple[int, int] | None:
    """Return the longest satellite-repeat block (centromere candidate), or None."""
    blocks = [b for b in kmer_blocks(sequence, satellite_motifs, max_gap=max_gap)
              if b[1] - b[0] >= min_array_len]
    if not blocks:
        return None
    return max(blocks, key=lambda b: b[1] - b[0])


def bin_chromosome(
    sequence: str,
    *,
    telomere_motifs: Sequence[str] = PLANT_TELOMERE_MOTIFS,
    satellite_motifs: Sequence[str] = (),
    end_window: int = 10_000,
    telomere_min_len: int = 100,
    centromere_min_len: int = 1000,
    pericentromere_flank: int = 50_000,
) -> list[StructureBin]:
    """Partition ``[0, len(sequence))`` into telomere/arm/pericentromere/centromere bins.

    Deterministic, threshold-based: capped ends become telomere bins, the longest
    satellite array becomes the centromere with flanking pericentromere, and the
    remainder is arm.  Returns contiguous, non-overlapping bins covering the sequence.
    (An optional learned arm/pericentromere classifier — sklearn — can refine this
    later; this rule-based binner needs no training and no extra dependency.)
    """
    n = len(sequence)
    if n == 0:
        return []

    boundaries: list[tuple[int, int, str]] = []

    telo = find_telomeres(
        sequence, telomere_motifs, end_window=end_window, min_array_len=telomere_min_len
    )
    left_telo = telo.five_prime_capped
    right_telo = telo.three_prime_capped
    telo_left_end = end_window if left_telo else 0
    telo_right_start = n - end_window if right_telo else n

    centro = find_centromere(
        sequence, satellite_motifs, min_array_len=centromere_min_len
    ) if satellite_motifs else None

    if left_telo:
        boundaries.append((0, telo_left_end, "telomere"))
    if right_telo:
        boundaries.append((telo_right_start, n, "telomere"))

    body_start = telo_left_end
    body_end = telo_right_start
    if centro is not None:
        c_start, c_end = centro
        c_start = max(body_start, c_start)
        c_end = min(body_end, c_end)
        peri_start = max(body_start, c_start - pericentromere_flank)
        peri_end = min(body_end, c_end + pericentromere_flank)
        if body_start < peri_start:
            boundaries.append((body_start, peri_start, "arm"))
        if peri_start < c_start:
            boundaries.append((peri_start, c_start, "pericentromere"))
        boundaries.append((c_start, c_end, "centromere"))
        if c_end < peri_end:
            boundaries.append((c_end, peri_end, "pericentromere"))
        if peri_end < body_end:
            boundaries.append((peri_end, body_end, "arm"))
    elif body_start < body_end:
        boundaries.append((body_start, body_end, "arm"))

    boundaries.sort()
    return [StructureBin(start=s, end=e, bin_type=t) for s, e, t in boundaries if e > s]
