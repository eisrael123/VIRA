"""Promoter-accessibility pipeline package."""

from typing import TYPE_CHECKING

__all__ = ["run_pipeline", "PipelineResult"]

if TYPE_CHECKING:
    from pipeline.run_pipeline import PipelineResult, run_pipeline


def __getattr__(name: str):
    if name in __all__:
        from pipeline import run_pipeline as _module

        return getattr(_module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
