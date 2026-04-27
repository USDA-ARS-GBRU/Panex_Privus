"""GFA adapter for the shared pangenome feature matrix."""

from __future__ import annotations

from privy.io.gfa import GfaGraph, get_samples_traversing_segment
from privy.pangenome.model import FeatureMatrix, FeatureRecord


def build_gfa_feature_matrix(graph: GfaGraph) -> FeatureMatrix:
    """Build a segment-level feature matrix from a parsed GFA graph."""
    all_samples = sorted(set(graph.sample_to_paths) | set(graph.sample_to_walks))
    features: list[FeatureRecord] = []
    presence: dict[str, frozenset[str]] = {}

    for seg_name in sorted(graph.segments):
        seg = graph.segments[seg_name]
        features.append(
            FeatureRecord(
                feature_id=seg.name,
                source_type="gfa",
                feature_type="segment",
                length=max(seg.length, 0),
                contig=seg.ref_contig,
                start=seg.ref_start,
                end=seg.ref_end,
            )
        )
        traversing = get_samples_traversing_segment(graph, seg_name)
        presence[seg.name] = frozenset(s for s in traversing if s in all_samples)

    return FeatureMatrix(
        source_type="gfa",
        features=tuple(features),
        samples=tuple(all_samples),
        presence=presence,
    )
