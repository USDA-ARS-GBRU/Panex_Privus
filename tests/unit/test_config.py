"""Unit tests for privy.core.config — YAML loading and validation."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from privy.core.config import (
    CohortConfig,
    PrivyConfig,
    ScanConfig,
    default_config,
    load_config,
)


class TestDefaultConfig:
    def test_returns_privy_config(self) -> None:
        cfg = default_config()
        assert isinstance(cfg, PrivyConfig)

    def test_default_mode(self) -> None:
        cfg = default_config()
        assert cfg.mode == "private_allele"

    def test_default_project_name(self) -> None:
        cfg = default_config()
        assert cfg.project_name == "privy_run"

    def test_default_min_target_support(self) -> None:
        cfg = default_config()
        assert cfg.scan.min_target_support == 1.0

    def test_default_max_off_target_support(self) -> None:
        cfg = default_config()
        assert cfg.scan.max_off_target_support == 0.0

    def test_scoring_defaults(self) -> None:
        cfg = default_config()
        assert cfg.scoring.discovery_weight == 1.0
        assert cfg.scoring.support_weight == 0.7
        assert cfg.scoring.penalty_weight == 0.8


class TestCohortConfigValidation:
    def test_overlap_raises(self) -> None:
        with pytest.raises(ValidationError):
            CohortConfig(targets=["S1", "S2"], off_targets=["S2", "S3"])

    def test_no_overlap_valid(self) -> None:
        cfg = CohortConfig(targets=["S1"], off_targets=["S2"])
        assert cfg.targets == ["S1"]


class TestLoadConfig:
    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_config(tmp_path / "nonexistent.yaml")

    def test_empty_yaml_returns_defaults(self, tmp_path: Path) -> None:
        yaml_file = tmp_path / "privy.yaml"
        yaml_file.write_text("")
        cfg = load_config(yaml_file)
        assert cfg.mode == "private_allele"

    def test_load_project_name(self, tmp_path: Path) -> None:
        yaml_file = tmp_path / "privy.yaml"
        yaml_file.write_text("project_name: my_soybean_run\n")
        cfg = load_config(yaml_file)
        assert cfg.project_name == "my_soybean_run"

    def test_load_cohort(self, tmp_path: Path) -> None:
        yaml_file = tmp_path / "privy.yaml"
        yaml_file.write_text(
            "cohorts:\n"
            "  targets: [Benning, Harosoy]\n"
            "  off_targets: [Jack, Lee]\n"
        )
        cfg = load_config(yaml_file)
        assert "Benning" in cfg.cohorts.targets
        assert "Jack" in cfg.cohorts.off_targets

    def test_load_scan_params(self, tmp_path: Path) -> None:
        yaml_file = tmp_path / "privy.yaml"
        yaml_file.write_text(
            "scan:\n"
            "  min_target_support: 0.8\n"
            "  max_off_target_support: 0.05\n"
            "  merge_distance: 500\n"
        )
        cfg = load_config(yaml_file)
        assert cfg.scan.min_target_support == pytest.approx(0.8)
        assert cfg.scan.merge_distance == 500

    def test_load_scoring_weights(self, tmp_path: Path) -> None:
        yaml_file = tmp_path / "privy.yaml"
        yaml_file.write_text(
            "scoring:\n"
            "  discovery_weight: 1.5\n"
            "  support_weight: 0.5\n"
        )
        cfg = load_config(yaml_file)
        assert cfg.scoring.discovery_weight == pytest.approx(1.5)

    def test_as_run_dict_serialisable(self, tmp_path: Path) -> None:
        import json

        cfg = default_config()
        d = cfg.as_run_dict()
        # Should not raise
        json.dumps(d)


class TestScanConfigValidation:
    def test_min_target_support_bounds(self) -> None:
        with pytest.raises(ValidationError):
            ScanConfig(min_target_support=1.5)
        with pytest.raises(ValidationError):
            ScanConfig(min_target_support=-0.1)

    def test_chunk_size_minimum(self) -> None:
        with pytest.raises(ValidationError):
            ScanConfig(chunk_size=100)  # below minimum of 1000
