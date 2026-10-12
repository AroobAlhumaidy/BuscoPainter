"""BuscoPainter: write BUSCO gene models onto an assembly as annotated GenBank.

The scientific logic lives in :mod:`buscopainter.painter` (internal version ``T7.1-fixed.3``):
the validated ``buscoPainter.T7.1.fixed.py`` (``T7.1-fixed.1``) plus the reporting and Artemis
changes of ``T7.1-fixed.2`` and the BUSCO-span correction of ``T7.1-fixed.3`` (CHANGELOG.md).
:mod:`buscopainter.validate` is the independent output validator; it does not import the painter.
"""

__version__ = "0.1.0"
REFERENCE_IMPLEMENTATION = "T7.1-fixed.3"
