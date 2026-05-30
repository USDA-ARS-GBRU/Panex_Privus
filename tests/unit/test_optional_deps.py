"""Unit tests for src/privy/utils/optional.py (tiered dependency policy)."""

from __future__ import annotations

import shutil

import pytest

from privy.utils.optional import (
    MissingDependencyError,
    MissingToolError,
    is_available,
    require,
    require_tool,
    tool_available,
)

# ---------------------------------------------------------------------------
# Tier 1: optional Python packages
# ---------------------------------------------------------------------------


class TestIsAvailable:
    def test_stdlib_module_available(self):
        assert is_available("json") is True

    def test_hard_dependency_available(self):
        assert is_available("numpy") is True   # core dependency

    def test_missing_module_not_available(self):
        assert is_available("definitely_not_a_real_module_xyz") is False


class TestRequire:
    def test_returns_module_when_present(self):
        mod = require("json", feature="x")
        assert mod.dumps({"a": 1}) == '{"a": 1}'

    def test_raises_with_feature_and_hint(self):
        with pytest.raises(MissingDependencyError) as exc:
            require("definitely_not_a_real_module_xyz", feature="DAPC structure")
        msg = str(exc.value)
        assert "DAPC structure" in msg
        assert "definitely_not_a_real_module_xyz" in msg
        assert "pip install" in msg

    def test_known_extra_uses_friendly_hint(self):
        # 'allel' maps to scikit-allel in the hint table; force the missing path
        # by requesting a guaranteed-absent name with a custom hint.
        with pytest.raises(MissingDependencyError) as exc:
            require("not_real_pkg", feature="f", hint="pip install special-thing")
        assert "special-thing" in str(exc.value)


# ---------------------------------------------------------------------------
# Tier 2: external CLI tools
# ---------------------------------------------------------------------------


class TestToolAvailable:
    def test_missing_tool_is_false(self):
        assert tool_available("definitely-not-a-real-tool-xyz") is False

    def test_present_tool_is_true(self):
        # Pick a tool that genuinely exists on this PATH; skip if none do.
        present = next((t for t in ("sh", "ls", "env", "python3") if shutil.which(t)), None)
        if present is None:
            pytest.skip("no common POSIX tool available to test against")
        assert tool_available(present) is True


class TestRequireTool:
    def test_raises_for_missing_tool(self):
        with pytest.raises(MissingToolError) as exc:
            require_tool("definitely-not-a-real-tool-xyz", feature="graph untangle")
        msg = str(exc.value)
        assert "graph untangle" in msg
        assert "definitely-not-a-real-tool-xyz" in msg
        assert "ingest its output" in msg

    def test_appends_hint(self):
        with pytest.raises(MissingToolError) as exc:
            require_tool("nope-xyz", feature="f", hint="conda install -c bioconda odgi")
        assert "conda install -c bioconda odgi" in str(exc.value)

    def test_returns_path_for_present_tool(self):
        present = next((t for t in ("sh", "ls", "env", "python3") if shutil.which(t)), None)
        if present is None:
            pytest.skip("no common POSIX tool available to test against")
        path = require_tool(present, feature="f")
        assert path == shutil.which(present)
