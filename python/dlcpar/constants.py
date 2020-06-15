"""
constants.py
Global constants for dlcpar
"""

#==========================================================
# MPR

DEFAULT_DUP_COST = 2
DEFAULT_LOSS_COST = 2
DEFAULT_COAL_COST = 1

#==========================================================
# landscapes

DEFAULT_DUP_RANGE = (0.2, 5)
DEFAULT_LOSS_RANGE = (0.2, 5)
DEFAULT_DUP_RANGE_STR = ':'.join(map(str, DEFAULT_DUP_RANGE))
DEFAULT_LOSS_RANGE_STR = ':'.join(map(str, DEFAULT_LOSS_RANGE))
