"""Small helpers for reading Fortran list-directed input and writing E/F edit descriptors
the way gfortran does, so legacy output files can be compared byte for byte."""
import math


class ListReader:
    """Mimics Fortran list-directed READ: each read starts a new record (line), consumes
    as many values as requested (continuing onto following lines if needed) and
    discards whatever is left on the last line it touched."""

    def __init__(self, path):
        with open(path) as f:
            self.lines = f.read().splitlines()
        self.pos = 0

    def line(self):
        s = self.lines[self.pos]
        self.pos += 1
        return s

    def read(self, n, types=None):
        vals = []
        while len(vals) < n:
            toks = self.line().replace(",", " ").split()
            vals.extend(toks)
        vals = vals[:n]
        if types is None:
            return vals
        if not isinstance(types, (list, tuple)):
            types = [types] * n
        return [_conv(v, t) for v, t in zip(vals, types)]


def _conv(tok, typ):
    if typ is int:
        return int(tok)
    return float(tok.replace("d", "e").replace("D", "e"))


def fe(x, w, d):
    """gfortran Ew.d."""
    if x == 0.0:
        mant, ex = "0." + "0" * d, 0
        neg = math.copysign(1.0, x) < 0
    elif math.isnan(x):
        return "NaN".rjust(w)
    else:
        neg = x < 0
        s = "%.*e" % (d - 1, abs(x))
        digits, e = s.split("e")
        digits = digits.replace(".", "")
        ex = int(e) + 1
        mant = "0." + digits
    if abs(ex) <= 99:
        es = "E%+03d" % ex
    else:
        es = "%+04d" % ex
    body = ("-" if neg else "") + mant + es
    if len(body) > w:
        body = body.replace("0.", ".", 1)
    if len(body) > w:
        return "*" * w
    return body.rjust(w)


def ff(x, w, d):
    """gfortran Fw.d."""
    s = "%.*f" % (d, x)
    if len(s) > w and (s.startswith("0.") or s.startswith("-0.")):
        s = s.replace("0.", ".", 1)
    if len(s) > w:
        return "*" * w
    return s.rjust(w)
