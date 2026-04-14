"""Configuration models and YAML loading for Panex Privus.

Config priority (lowest → highest):
    1. Package defaults (field defaults in Pydantic models)
    2. YAML config file (``--config`` flag)
    3. CLI overrides (applied by each subcommand after loading)

The resolved configuration must be written to ``run.json`` so every run is
fully reproducible.

Example usage::

    from privy.core.config import load_config, default_config

    cfg = load_config(Path("privy.yaml"))
    # CLI override:
    cfg = cfg.model_copy(update={"scan": cfg.scan.model_copy(update={"merge_distance": 500})})
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal, Optional

import yaml
from pydantic import BaseModel, Field, model_validator


class CohortConfig(BaseModel):
    """Cohort sample assignments from YAML config."""

    targets: list[str] = Field(default_factory=list, description="Target sample names.")
    off_targets: list[str] = Field(
        default_factory=list, description="Off-target sample names."
    )
    ignored_samples: list[str] = Field(
        default_factory=list, description="Samples to exclude from analysis."
    )

    @model_validator(mode="after")
    def check_no_overlap(self) -> "CohortConfig":
        overlap = set(self.targets) & set(self.off_targets)
        if overlap:
            raise ValueError(
                f"Samples appear in both targets and off_targets: {sorted(overlap)}"
            )
        return self


class ScanConfig(BaseModel):
    """Discovery parameters for ``privy scan``."""

    min_target_support: float = Field(
        default=1.0, ge=0.0, le=1.0,
        description="Minimum fraction of target samples that must support the allele.",
    )
    max_off_target_support: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Maximum fraction of off-target samples allowed to support the allele.",
    )
    merge_distance: int = Field(
        default=0, ge=0,
        description="Merge loci within this bp distance into candidate regions. 0 = no merge.",
    )
    strictness_report: bool = Field(
        default=True,
        description="Include strictness_class in all locus outputs.",
    )
    relaxed_target_missing: Optional[float] = Field(
        default=None, ge=0.0, le=1.0,
        description="Tolerated target missingness fraction for relaxed_threshold class.",
    )
    relaxed_offtarget_missing: Optional[float] = Field(
        default=None, ge=0.0, le=1.0,
        description="Tolerated off-target missingness fraction for relaxed_threshold class.",
    )
    chunk_size: int = Field(
        default=100_000, ge=1000,
        description="Contig chunk size (bp) for streaming large inputs.",
    )
    same_variant_class_only: bool = Field(
        default=False,
        description="Only merge loci of the same variant class into regions.",
    )
    pass_only: bool = Field(
        default=True,
        description="Require VCF FILTER=PASS.",
    )
    min_qual: Optional[float] = Field(
        default=None,
        description="Minimum VCF QUAL score.",
    )
    allow_multiallelic: bool = Field(
        default=True,
        description="Evaluate multiallelic VCF records.",
    )
    mode: Literal["private_allele", "private_genotype", "private_sv_state"] = Field(
        default="private_allele",
        description="Discovery mode.",
    )


class BamConfig(BaseModel):
    """Parameters for BAM support layer."""

    enabled: bool = False
    min_depth: int = Field(default=8, ge=0)
    allele_fraction_min: float = Field(default=0.2, ge=0.0, le=1.0)
    min_alt_count: int = Field(default=2, ge=0)
    summarize_softclips: bool = False
    summarize_splitreads: bool = False


class GfaConfig(BaseModel):
    """Parameters for GFA graph-context and standalone GFA scan."""

    enabled: bool = False
    junction_window_bp: int = Field(default=1000, ge=0)
    report_path_membership: bool = True
    report_graph_complexity: bool = True
    min_segment_length: int = Field(
        default=1, ge=1,
        description=(
            "Minimum segment length (bp) to evaluate.  Shorter segments are "
            "skipped.  Increase to filter out very short bubbles."
        ),
    )
    path_name_format: str = Field(
        default="pangenome",
        description=(
            "Convention for extracting sample names from P-line path names. "
            "'pangenome' = SAMPLE#HAP#CONTIG (minigraph-cactus / PGGB); "
            "'plain' = the full path name is the sample name."
        ),
    )


class XmfaConfig(BaseModel):
    """Parameters for XMFA alignment-corroboration layer."""

    enabled: bool = False
    gap_aware: bool = True
    window_bp: int = Field(default=500, ge=0)


class CompareConfig(BaseModel):
    """Parameters for ``privy compare``."""

    overlap_mode: Literal["any", "reciprocal", "contained"] = "reciprocal"
    min_reciprocal_overlap: float = Field(default=0.5, ge=0.0, le=1.0)
    breakpoint_tolerance_bp: int = Field(default=200, ge=0)
    require_state_compatibility: bool = False


class ScoringConfig(BaseModel):
    """Scoring weights.  Stored in run.json for reproducibility."""

    discovery_weight: float = Field(default=1.0, ge=0.0)
    support_weight: float = Field(default=0.7, ge=0.0)
    penalty_weight: float = Field(default=0.8, ge=0.0)


class PrivyConfig(BaseModel):
    """Top-level configuration model for a Panex Privus run.

    All fields have sensible defaults so the package works with zero
    configuration.  A YAML file and CLI flags each override the previous
    layer.
    """

    project_name: str = "privy_run"
    mode: Literal["private_allele", "private_genotype", "private_sv_state"] = "private_allele"
    cohorts: CohortConfig = Field(default_factory=CohortConfig)
    scan: ScanConfig = Field(default_factory=ScanConfig)
    bam: BamConfig = Field(default_factory=BamConfig)
    gfa: GfaConfig = Field(default_factory=GfaConfig)
    xmfa: XmfaConfig = Field(default_factory=XmfaConfig)
    compare: CompareConfig = Field(default_factory=CompareConfig)
    scoring: ScoringConfig = Field(default_factory=ScoringConfig)

    def as_run_dict(self) -> dict:  # type: ignore[type-arg]
        """Return a JSON-serialisable dict for writing into run.json."""
        return self.model_dump(mode="json")


def load_config(path: Path) -> PrivyConfig:
    """Load a YAML config file and return a validated :class:`PrivyConfig`.

    Args:
        path: Path to a ``privy.yaml`` file.

    Returns:
        A fully validated :class:`PrivyConfig` instance.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        pydantic.ValidationError: If the YAML contains invalid values.
    """
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    if raw is None:
        raw = {}
    return PrivyConfig.model_validate(raw)


def default_config() -> PrivyConfig:
    """Return a :class:`PrivyConfig` with all package defaults."""
    return PrivyConfig()
