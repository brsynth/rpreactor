from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# Expose main classes for end-user convenience
from rpreactor.rule.burner import RuleBurner
