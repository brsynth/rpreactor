from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# Expose main classes for end-user convenience
from rpreactor.rule.burner import RuleBurner
from rpreactor.rule.utils import create_db_from_retrorules


# Silence most of RDKit complaints
from rdkit import RDLogger

RD_LOGGER = RDLogger.logger()
RD_LOGGER.setLevel(RDLogger.CRITICAL)
