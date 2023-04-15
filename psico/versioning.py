import re


def make_version_int_tuple(version: str) -> tuple:
    """
    Parse a version string into an integer tuple without trailing zeros.

    Ignores all additional qualifiers like "alpha" or "rc".
    """
    m = re.match(r"[0-9.]*", version)
    parts = [int(p) for p in m.group().split(".") if p]
    while parts and parts[-1] == 0:
        parts.pop()
    return tuple(parts)
