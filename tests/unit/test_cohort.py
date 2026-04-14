"""Unit tests for privy.core.cohort.CohortDefinition."""

import pytest

from privy.core.cohort import CohortDefinition


class TestCohortDefinitionConstruction:
    def test_basic_construction(self) -> None:
        cohort = CohortDefinition.from_lists(
            targets=["S1", "S2"],
            off_targets=["S3", "S4"],
        )
        assert cohort.n_targets == 2
        assert cohort.n_off_targets == 2
        assert cohort.ignored_samples == ()

    def test_with_ignored_samples(self) -> None:
        cohort = CohortDefinition.from_lists(
            targets=["S1"],
            off_targets=["S2"],
            ignored_samples=["S3", "S4"],
        )
        assert len(cohort.ignored_samples) == 2

    def test_with_metadata(self) -> None:
        cohort = CohortDefinition.from_lists(
            targets=["S1"],
            off_targets=["S2"],
            metadata={"species": "Glycine max"},
        )
        assert cohort.metadata["species"] == "Glycine max"

    def test_direct_construction(self) -> None:
        cohort = CohortDefinition(
            targets=("A", "B"),
            off_targets=("C",),
        )
        assert cohort.n_targets == 2


class TestCohortDefinitionValidation:
    def test_overlap_raises(self) -> None:
        with pytest.raises(ValueError, match="both targets and off_targets"):
            CohortDefinition.from_lists(
                targets=["S1", "S2"],
                off_targets=["S2", "S3"],
            )

    def test_empty_targets_raises(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            CohortDefinition.from_lists(targets=[], off_targets=["S1"])

    def test_empty_off_targets_raises(self) -> None:
        with pytest.raises(ValueError, match="at least one"):
            CohortDefinition.from_lists(targets=["S1"], off_targets=[])

    def test_multi_overlap_error_lists_all(self) -> None:
        with pytest.raises(ValueError) as exc_info:
            CohortDefinition.from_lists(
                targets=["A", "B", "C"],
                off_targets=["B", "C", "D"],
            )
        msg = str(exc_info.value)
        assert "B" in msg
        assert "C" in msg


class TestCohortDefinitionLookups:
    @pytest.fixture
    def cohort(self) -> CohortDefinition:
        return CohortDefinition.from_lists(
            targets=["T1", "T2"],
            off_targets=["O1", "O2"],
            ignored_samples=["X1"],
        )

    def test_is_target(self, cohort: CohortDefinition) -> None:
        assert cohort.is_target("T1")
        assert cohort.is_target("T2")
        assert not cohort.is_target("O1")

    def test_is_off_target(self, cohort: CohortDefinition) -> None:
        assert cohort.is_off_target("O1")
        assert not cohort.is_off_target("T1")

    def test_is_ignored(self, cohort: CohortDefinition) -> None:
        assert cohort.is_ignored("X1")
        assert not cohort.is_ignored("T1")

    def test_is_known(self, cohort: CohortDefinition) -> None:
        assert cohort.is_known("T1")
        assert cohort.is_known("O2")
        assert cohort.is_known("X1")
        assert not cohort.is_known("UNKNOWN")

    def test_all_active_samples(self, cohort: CohortDefinition) -> None:
        active = cohort.all_active_samples
        assert active == {"T1", "T2", "O1", "O2"}
        assert "X1" not in active

    def test_repr(self, cohort: CohortDefinition) -> None:
        r = repr(cohort)
        assert "n_targets=2" in r
        assert "n_off_targets=2" in r


class TestCohortDefinitionImmutability:
    def test_frozen(self) -> None:
        cohort = CohortDefinition.from_lists(["T1"], ["O1"])
        with pytest.raises((AttributeError, TypeError)):
            cohort.targets = ("T2",)  # type: ignore[misc]
