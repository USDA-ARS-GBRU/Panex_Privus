"""Deterministic synthetic pangenome generator for tests, CI, and examples.

Builds tiny (KB-scale) GFA graphs that *look like* real crop pangenome output —
PanSN-named paths (``sample#hap#contig``), rGFA ``SN``/``SO`` tags on the
reference's segments, target/off-target cohort labels, and optionally planted
structural events (inversion, duplication, translocation) — so the synteny,
projection, microhaplotype, and population-genetics layers can be developed and
tested on realistic-shaped inputs *without downloading anything*.

Everything is deterministic: identical inputs always produce byte-identical GFA,
so fixtures are reproducible and diff-stable.  Pure-Python, no dependencies.

Real-world (GB-scale) validation happens separately on the UGA Sapelo2 cluster;
see the development roadmap.

Coordinate convention: 0-based, half-open — matching the rest of Privy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

_BASES = "ACGT"


def _deterministic_sequence(seg_id: str, length: int) -> str:
    """A stable, segment-specific DNA sequence of the requested length."""
    # Seed the rotation from the segment id so different segments differ, but the
    # output is fully deterministic (no RNG).
    offset = sum(ord(c) for c in seg_id) % len(_BASES)
    rotated = _BASES[offset:] + _BASES[:offset]
    repeats = (rotated * (length // len(rotated) + 1))[:length]
    return repeats


@dataclass
class _Segment:
    seg_id: str
    sequence: str
    ref_contig: str | None = None   # SN:Z (only on the reference's segments)
    ref_start: int | None = None    # SO:i


@dataclass
class _Genome:
    name: str                       # full PanSN path name sample#hap#contig
    steps: list[tuple[str, str]]    # (segment_id, orientation)
    cohort: str                     # "target" or "offtarget"


@dataclass
class SyntheticPangenome:
    """Builder for a small synthetic GFA pangenome.

    Add segments, then genomes (each a walk of ``(segment, orientation)`` steps),
    then render to GFA text or write to a file.
    """

    segments: dict[str, _Segment] = field(default_factory=dict)
    genomes: list[_Genome] = field(default_factory=list)

    # -- construction ------------------------------------------------------

    def add_segment(self, seg_id: str, length: int) -> SyntheticPangenome:
        """Add a segment with a deterministic sequence of *length* bp."""
        if seg_id in self.segments:
            raise ValueError(f"duplicate segment id {seg_id!r}")
        self.segments[seg_id] = _Segment(seg_id, _deterministic_sequence(seg_id, length))
        return self

    def add_genome(
        self,
        name: str,
        steps: list[tuple[str, str]],
        *,
        cohort: str = "offtarget",
    ) -> SyntheticPangenome:
        """Add a genome as an ordered walk of ``(segment_id, orientation)`` steps."""
        for seg_id, orient in steps:
            if seg_id not in self.segments:
                raise ValueError(f"genome {name!r} references unknown segment {seg_id!r}")
            if orient not in ("+", "-"):
                raise ValueError(f"orientation must be '+' or '-', got {orient!r}")
        self.genomes.append(_Genome(name=name, steps=list(steps), cohort=cohort))
        return self

    def tag_reference(self, genome_name: str) -> SyntheticPangenome:
        """Assign rGFA ``SN``/``SO`` tags to segments along *genome_name*'s path.

        Mimics minigraph's stable-coordinate model: segments first seen on the
        reference get a stable ``(contig, offset)``; others stay untagged.
        """
        genome = next((g for g in self.genomes if g.name == genome_name), None)
        if genome is None:
            raise ValueError(f"unknown genome {genome_name!r}")
        contig = genome.name.rsplit("#", 1)[-1]
        cursor = 0
        for seg_id, _orient in genome.steps:
            seg = self.segments[seg_id]
            if seg.ref_contig is None:   # first occurrence only
                seg.ref_contig = contig
                seg.ref_start = cursor
            cursor += len(seg.sequence)
        return self

    # -- cohort accessors --------------------------------------------------

    def cohort(self, which: str) -> list[str]:
        """Return genome names in cohort *which* ("target" / "offtarget")."""
        return [g.name for g in self.genomes if g.cohort == which]

    # -- rendering ---------------------------------------------------------

    def to_gfa(self, *, use_walks: bool = False) -> str:
        """Render the pangenome to GFA1.1 text.

        Args:
            use_walks: Emit genomes as W-lines (walks) instead of P-lines (paths).
        """
        lines: list[str] = ["H\tVN:Z:1.1"]

        for seg in self.segments.values():
            fields = [f"S\t{seg.seg_id}\t{seg.sequence}\tLN:i:{len(seg.sequence)}"]
            if seg.ref_contig is not None:
                fields.append(f"SN:Z:{seg.ref_contig}")
                fields.append(f"SO:i:{seg.ref_start}")
            lines.append("\t".join(fields))

        # Deduplicated links between consecutive steps (orientation-aware).
        seen_links: set[tuple[str, str, str, str]] = set()
        link_lines: list[str] = []
        for genome in self.genomes:
            for (a_seg, a_orient), (b_seg, b_orient) in zip(
                genome.steps, genome.steps[1:], strict=False
            ):
                key = (a_seg, a_orient, b_seg, b_orient)
                if key in seen_links:
                    continue
                seen_links.add(key)
                link_lines.append(f"L\t{a_seg}\t{a_orient}\t{b_seg}\t{b_orient}\t0M")
        lines.extend(link_lines)

        for genome in self.genomes:
            if use_walks:
                lines.append(self._walk_line(genome))
            else:
                lines.append(self._path_line(genome))

        return "\n".join(lines) + "\n"

    def _path_line(self, genome: _Genome) -> str:
        seg_field = ",".join(f"{seg}{orient}" for seg, orient in genome.steps)
        return f"P\t{genome.name}\t{seg_field}\t*"

    def _walk_line(self, genome: _Genome) -> str:
        sample, hap, contig = _split_name(genome.name)
        total = sum(len(self.segments[seg].sequence) for seg, _ in genome.steps)
        walk = "".join(
            f"{'>' if orient == '+' else '<'}{seg}" for seg, orient in genome.steps
        )
        return f"W\t{sample}\t{hap}\t{contig}\t0\t{total}\t{walk}"

    def write(self, path: Path, *, use_walks: bool = False) -> Path:
        """Write the GFA to *path* and return it."""
        path = Path(path)
        path.write_text(self.to_gfa(use_walks=use_walks), encoding="utf-8")
        return path


def _split_name(name: str) -> tuple[str, int, str]:
    parts = name.split("#")
    if len(parts) >= 3:
        try:
            hap = int(parts[1])
        except ValueError:
            hap = 0
        return parts[0], hap, parts[-1]
    return name, 0, name


# ---------------------------------------------------------------------------
# Canonical scenario builders
# ---------------------------------------------------------------------------


def collinear_pangenome(
    n_genomes: int = 4,
    n_segments: int = 6,
    seg_len: int = 10,
    *,
    n_target: int = 2,
    contig: str = "chr1",
) -> SyntheticPangenome:
    """All genomes traverse the same segments in the same order (perfect collinearity).

    The first *n_target* genomes are labelled ``target``, the rest ``offtarget``.
    The first genome is tagged as the rGFA reference.
    """
    pg = SyntheticPangenome()
    seg_ids = [f"s{i}" for i in range(1, n_segments + 1)]
    for seg_id in seg_ids:
        pg.add_segment(seg_id, seg_len)
    steps = [(seg_id, "+") for seg_id in seg_ids]
    for g in range(n_genomes):
        cohort = "target" if g < n_target else "offtarget"
        pg.add_genome(f"sample{g}#0#{contig}", steps, cohort=cohort)
    pg.tag_reference(f"sample0#0#{contig}")
    return pg


def inversion_pangenome(seg_len: int = 10, *, contig: str = "chr1") -> SyntheticPangenome:
    """4 collinear genomes; the last genome carries an inverted middle run (s3,s4)."""
    pg = SyntheticPangenome()
    seg_ids = [f"s{i}" for i in range(1, 7)]
    for seg_id in seg_ids:
        pg.add_segment(seg_id, seg_len)
    forward = [(s, "+") for s in seg_ids]
    # invert s3,s4: reverse their order and flip orientation
    inverted = (
        [("s1", "+"), ("s2", "+"), ("s4", "-"), ("s3", "-"), ("s5", "+"), ("s6", "+")]
    )
    for g in range(3):
        pg.add_genome(f"sample{g}#0#{contig}", forward, cohort="offtarget")
    pg.add_genome(f"sample3#0#{contig}", inverted, cohort="target")
    pg.tag_reference(f"sample0#0#{contig}")
    return pg


def duplication_pangenome(seg_len: int = 10, *, contig: str = "chr1") -> SyntheticPangenome:
    """3 collinear genomes; one target genome has a tandem duplication of s2."""
    pg = SyntheticPangenome()
    seg_ids = [f"s{i}" for i in range(1, 5)]
    for seg_id in seg_ids:
        pg.add_segment(seg_id, seg_len)
    forward = [(s, "+") for s in seg_ids]
    duplicated = [("s1", "+"), ("s2", "+"), ("s2", "+"), ("s3", "+"), ("s4", "+")]
    pg.add_genome(f"sample0#0#{contig}", forward, cohort="offtarget")
    pg.add_genome(f"sample1#0#{contig}", forward, cohort="offtarget")
    pg.add_genome(f"sample2#0#{contig}", duplicated, cohort="target")
    pg.tag_reference(f"sample0#0#{contig}")
    return pg


def allopolyploid_pangenome(seg_len: int = 10) -> SyntheticPangenome:
    """An AADD-style allotetraploid: two subgenomes (chrA, chrD) per sample.

    Each sample contributes two haplotype paths — ``sample#0#chrA`` and
    ``sample#1#chrD`` — over disjoint segment sets, so subgenome phasing and
    homeolog-aware logic can be exercised.
    """
    pg = SyntheticPangenome()
    a_segs = [f"a{i}" for i in range(1, 5)]
    d_segs = [f"d{i}" for i in range(1, 5)]
    for seg_id in a_segs + d_segs:
        pg.add_segment(seg_id, seg_len)
    a_steps = [(s, "+") for s in a_segs]
    d_steps = [(s, "+") for s in d_segs]
    for g in range(3):
        cohort = "target" if g == 0 else "offtarget"
        pg.add_genome(f"sample{g}#0#chrA", a_steps, cohort=cohort)
        pg.add_genome(f"sample{g}#1#chrD", d_steps, cohort=cohort)
    pg.tag_reference("sample0#0#chrA")
    pg.tag_reference("sample0#1#chrD")
    return pg
