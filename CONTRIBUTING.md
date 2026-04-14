# Contributing to Panex Privus

Thank you for your interest in contributing. This document covers how to report issues, set up a development environment, and submit changes.

---

## Reporting bugs and requesting features

Please use [GitHub Issues](https://github.com/USDA-ARS-GBRU/Panex_Privus/issues).

When reporting a bug, include:

- The exact command you ran (redact any private sample names if needed)
- The full error message or traceback
- Your operating system and Python version (`python --version`)
- The Panex Privus version (`privy --version`)
- A small example VCF if possible — the smaller the better

When requesting a feature, describe:

- The biological or analytical problem you are trying to solve
- What inputs you have and what output you expect
- Whether this is blocking your research or a nice-to-have

---

## Development setup

```bash
git clone https://github.com/USDA-ARS-GBRU/Panex_Privus.git
cd Panex_Privus

# Create and activate a virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate   # on Windows: .venv\Scripts\activate

# Install in editable mode with development dependencies
pip install -e ".[dev]"
```

Verify the setup:

```bash
pytest             # run all tests
ruff check src/    # lint
mypy src/privy/core/   # type-check the core module
```

---

## Running the tests

```bash
# All tests
pytest

# Unit tests only (fast, no pysam required for most)
pytest tests/unit/

# Integration tests (requires pysam)
pytest tests/integration/

# With coverage report
pytest --cov=privy --cov-report=term-missing
```

All 183 tests must pass before a pull request is merged. New features require new tests.

---

## Code style

This project uses:

- **[ruff](https://docs.astral.sh/ruff/)** for linting and import sorting
- **[mypy](https://mypy.readthedocs.io/)** for static type checking (core module)
- Standard Python docstring style (Google-style within each function/class)

Run before committing:

```bash
ruff check src/ tests/
ruff format src/ tests/
mypy src/privy/core/
```

CI enforces these checks on every pull request.

---

## Project architecture

Before making changes, read [`docs/architecture.md`](docs/architecture.md). The key principle:

> File formats do not define truth. Each source (VCF, BAM, GFA, XMFA) contributes evidence into a common internal representation. The core logic (`src/privy/core/`) must remain format-independent.

The `StrictnessClass` framework is a deliberate design choice — missingness must never be silently folded into pass/fail. Do not simplify it.

### Directory layout

```
src/privy/
├── cli/          — typer subcommands; argument parsing only
├── core/         — domain objects and pure logic (no file I/O)
├── io/           — format readers and writers
├── backends/     — format-specific scan orchestrators
├── compare/      — cross-evidence comparison logic
├── report/       — report generation (stub)
├── plot/         — visualization (stub)
└── utils/        — logging, config, metrics, misc
```

---

## Submitting a pull request

1. Fork the repository and create a branch from `main`
2. Make your changes
3. Add or update tests for any changed behavior
4. Run the full test suite (`pytest`) and linters (`ruff`, `mypy`)
5. Write a clear PR description explaining what changed and why
6. Open the pull request

For large architectural changes, please open an issue for discussion first — it saves everyone time.

---

## Development roadmap

See the [roadmap in README.md](README.md#current-status) for the planned phases. If you want to contribute to a specific phase, comment on the relevant issue or open a new one.

Areas currently most in need of contribution:

- BAM support layer (`src/privy/backends/bam_support.py`) — Phase 3
- Report generation (`src/privy/report/`) — Phase 2
- Plot command (`src/privy/plot/`) — Phase 2
- Additional integration test fixtures

---

## Code of conduct

Be respectful and constructive. This project is affiliated with academic research; all contributions should meet the standards expected in a professional scientific environment.
