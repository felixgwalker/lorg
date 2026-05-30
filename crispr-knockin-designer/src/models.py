"""Data models for CRISPR Knock-in Designer."""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum


class DonorType(Enum):
    SSODDN = "ssODN"
    DSDNA = "dsDNA"
    AAV = "AAV"


@dataclass
class HDRTemplate:
    left_arm: str
    insert_sequence: str
    right_arm: str
    pam_mutation: str
    donor_sequence: str
    donor_type: DonorType
    total_length: int


@dataclass
class KnockinDesign:
    guide_spacer: str
    pam: str
    cut_position: int
    hdr_template: HDRTemplate
    output_files: dict[str, str] = field(default_factory=dict)
    pipeline_version: str = "1.0.0"
