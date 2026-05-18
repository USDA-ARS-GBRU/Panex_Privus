"""Interactive self-contained HTML dashboards for Panex Privus."""

from privy.interactive.focus import run_focus_dashboards
from privy.interactive.models import FocusRegion, parse_focus_region

__all__ = ["FocusRegion", "parse_focus_region", "run_focus_dashboards"]
