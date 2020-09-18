"""
Exceptions raised when something went wrong processing a rule.
"""


class ChemConversionError(Exception):
    """Raised when something went wrong during chemical conversion to RDKit mol object."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "CHEM-CONVERSION-ERROR: {}".format(self._msg)


class RuleConversionError(Exception):
    """Raised when something went wrong during SMARTS conversion to RDKit rxn object."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-CONVERSION-ERROR: {}".format(self._msg)


class RuleMatchError(Exception):
    """Raised when something went wrong when matching a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-MATCH-ERROR: {}".format(self._msg)


class RuleFireError(Exception):
    """Raised when something went wrong when firing a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-FIRE-ERROR: {}".format(self._msg)
