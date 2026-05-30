"""Data models for CRISPR Delivery Strategy Selector."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class IntegrationRisk(Enum):
    NONE = "none"
    VERY_LOW = "very_low"
    LOW = "low"
    HIGH = "high"


@dataclass
class DeliveryStrategy:
    name: str
    description: str
    compatibility_score: float
    integration_risk: IntegrationRisk
    payload_size_ok: bool
    notes: list[str] = field(default_factory=list)


@dataclass
class DeliveryScore:
    strategies: list[DeliveryStrategy]
    recommended: str
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
