from __future__ import annotations


class MartiniMapperError(Exception):
    """Base exception for martini_mapper."""


class InvalidSmilesError(MartiniMapperError, ValueError):
    """Raised when a SMILES string is empty or cannot be parsed by RDKit."""


class MappingError(MartiniMapperError):
    """Raised when the mapping algorithm cannot assign beads consistently."""


class BeadNotFoundError(MappingError, KeyError):
    """Raised when the dictionary contains no bead that matches a required fragment."""


class BeadAssignmentError(MappingError):
    """Raised when a bead assignment would be inconsistent (e.g., overlaps)."""


class ExternalToolError(MartiniMapperError):
    """Raised when an external dependency/tooling step fails (e.g., xtb)."""


class OutputError(MartiniMapperError):
    """Raised when writing or organizing output files fails."""


class ConfigError(MartiniMapperError):
    """Raised when configuration is invalid or missing."""
