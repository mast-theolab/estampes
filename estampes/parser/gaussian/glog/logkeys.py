"""Store patterns for common keywords in Gaussian log file.

To limit redundancy, common or complex patters are stored as variables.
"""

# RR_OMEGA = r'^\s+-> (?P<val>Omega =\s*\d+\.\d+ cm.-1.*?)\s*$'
# RR_OMEGA_LINE = r'^\s+-> (?P<val>Omega =\s*\d+\.\d+ cm.-1' \
#     + r'(?:, Gamma =\s*\d+\.\d+ cm.-1))?\s*$'
RR_OMEGA_LINE = r'^\s+-> (?P<val>Omega =\s*\d+\.\d+ cm.-1.*?)\s*$'
RR_OMEGA_UNIT = r'^\s+-> Omega = \s*(?P<val>\d+\.\d+ cm.-1)' \
    + r'(, Gamma =\s*\d+\.\d+ cm.-1)?\s*$'
RR_OMEGA_VAL = r'^\s+-> Omega = \s*(?P<val>\d+\.\d+) cm.-1' \
    + r'(, Gamma =\s*\d+\.\d+ cm.-1)?\s*$'
