from .peak import PeakFrame
from collections import Counter
import numpy as np
import re
from typing import Dict


EXACT_MASS = {
    "e":5.4858E-4,
    # 1–10
    "H": 1.00782503223,      # 1H
    "He": 4.00260325413,     # 4He
    "Li": 7.0160034366,      # 7Li
    "Be": 9.012183065,       # 9Be
    "B": 11.00930536,        # 11B
    "C": 12.00000000000,     # 12C (定义)
    "N": 14.00307400443,     # 14N
    "O": 15.99491461957,     # 16O
    "F": 18.99840316273,     # 19F
    "Ne": 19.9924401762,     # 20Ne

    # 11–20
    "Na": 22.9897692820,     # 23Na
    "Mg": 23.985041697,      # 24Mg
    "Al": 26.98153853,       # 27Al
    "Si": 27.97692653465,    # 28Si
    "P": 30.97376199842,     # 31P
    "S": 31.9720711744,      # 32S
    "Cl": 34.968852682,      # 35Cl
    "Ar": 35.967545105,      # 36Ar
    "K": 38.9637064864,      # 39K
    "Ca": 39.962590863,      # 40Ca

    # 21–30
    "Sc": 44.95590828,       # 45Sc
    "Ti": 47.94794198,       # 48Ti
    "V": 50.94395704,        # 51V
    "Cr": 51.94050623,       # 52Cr
    "Mn": 54.93804391,       # 55Mn
    "Fe": 55.93493633,       # 56Fe
    "Co": 58.93319429,       # 59Co
    "Ni": 57.93534241,       # 58Ni (最丰之一)
    "Cu": 62.92959772,       # 63Cu
    "Zn": 63.92914201,       # 64Zn (较丰)

    # 31–40
    "Ga": 68.9255735,        # 69Ga
    "Ge": 73.921177761,      # 74Ge
    "As": 74.92159457,       # 75As
    "Se": 79.9165218,        # 80Se (丰度较高)
    "Br": 78.9183376,        # 79Br
    "Kr": 83.9114977282,     # 84Kr
    "Rb": 84.9117897379,     # 85Rb
    "Sr": 87.9056125,        # 88Sr
    "Y": 88.9058403,         # 89Y
    "Zr": 89.9046977,        # 90Zr

    # 41–50
    "Nb": 92.9063730,        # 93Nb
    "Mo": 97.90540482,       # 98Mo
    "Tc": 98.9062547,        # 99Tc (长寿命放射性)
    "Ru": 101.9043441,       # 102Ru
    "Rh": 102.9054980,       # 103Rh
    "Pd": 105.9034804,       # 106Pd
    "Ag": 106.9050916,       # 107Ag
    "Cd": 113.9033651,       # 114Cd
    "In": 114.9038788,       # 115In
    "Sn": 119.90220163,      # 120Sn

    # 51–60
    "Sb": 120.9038120,       # 121Sb
    "Te": 129.9062244,       # 130Te
    "I": 126.9044719,        # 127I
    "Xe": 131.90415509,      # 132Xe
    "Cs": 132.905451961,     # 133Cs
    "Ba": 137.9052472,       # 138Ba
    "La": 138.9063563,       # 139La
    "Ce": 139.9054431,       # 140Ce
    "Pr": 140.9076576,       # 141Pr
    "Nd": 141.9077290,       # 142Nd

    # 61–70
    "Pm": 144.9127559,       # 145Pm (放射性)
    "Sm": 151.9197397,       # 152Sm
    "Eu": 152.9212380,       # 153Eu
    "Gd": 157.9241123,       # 158Gd
    "Tb": 158.9253547,       # 159Tb
    "Dy": 163.9291819,       # 164Dy
    "Ho": 164.9303288,       # 165Ho
    "Er": 165.9302995,       # 166Er
    "Tm": 168.9342179,       # 169Tm
    "Yb": 173.9388664,       # 174Yb

    # 71–80
    "Lu": 174.9407752,       # 175Lu
    "Hf": 179.9465570,       # 180Hf
    "Ta": 180.9479958,       # 181Ta
    "W": 183.9509309,        # 184W
    "Re": 186.9557501,       # 187Re
    "Os": 191.9614770,       # 192Os
    "Ir": 192.9629216,       # 193Ir
    "Pt": 194.9647917,       # 195Pt
    "Au": 196.9665688,       # 197Au
    "Hg": 201.9706434,       # 202Hg

    # 81–90
    "Tl": 204.9744278,       # 205Tl
    "Pb": 207.9766525,       # 208Pb
    "Bi": 208.9803991,       # 209Bi
    "Po": 208.9824308,       # 209Po (放射性)
    "At": 209.9871479,       # 210At (放射性)
    "Rn": 222.0175777,       # 222Rn (放射性)
    "Fr": 223.0197359,       # 223Fr (放射性)
    "Ra": 226.0254098,       # 226Ra (放射性)
    "Ac": 227.0277523,       # 227Ac
    "Th": 232.0380558,       # 232Th

    # 91–100
    "Pa": 231.0358842,       # 231Pa
    "U": 238.0507884,        # 238U
    "Np": 237.0481736,       # 237Np
    "Pu": 244.064204,        # 244Pu
    "Am": 243.0613813,       # 243Am
    "Cm": 247.070354,        # 247Cm
    "Bk": 247.070307,        # 247Bk
    "Cf": 251.0795886,       # 252Cf (近似)
    "Es": 252.08298,         # 252Es (近似)
    "Fm": 257.0951061,       # 257Fm (近似)

    # 101–108
    "Md": 258.0984315,       # 258Md (近似)
    "No": 259.10103,         # 259No (近似)
    "Lr": 262.10961,         # 262Lr (近似)
    "Rf": 267.12179,         # 267Rf (近似)
    "Db": 268.12567,         # 268Db (近似)
    "Sg": 271.13393,         # 271Sg (近似)
    "Bh": 272.13826,         # 272Bh (近似)
    "Hs": 270.13429,         # 270Hs (近似)
}

## C2H6OS = DMSO    C3H8O = iPrOH

ION_POS = [
    "[M+H-C12H20O9]+",  "[M-C16H17N2O3]+",      "[M+H-C9H10O5]+",       "[M-C15H25+Na]+",       "[M-C6H10ClO4]+",
    "[M-C6H11O6]+",      "[M-C9H5O4]+",          "[M-C6H8O6+H]+",        "[M-C6H8NO5]+",         "[M-C14H23+Na]+",
    "[M-C6H10O5+H]+",    "[M-C5H8NO4]+",         "[M-C9H7O2+H]+",        "[M-C6H10O4+H]+",       "[M-C8H10NO]+",
    "[M-C5H4N5]+",       "[M-C11H21+Na]+",       "[M-C8H17O]+",          "[M-C6H6N2O+H]+",       "[M-C10H17+Na]+",
    "[M-C7H13O]+",       "[M-C4H3N2O2]+",        "[M-C9H15+Na]+",        "[M-C7H15]+",           "[M-C6H6N2O+Na]+",
    "[M-H2O4P]+",        "[M-5H2O+H]+",          "[M-C3H7N2O]+",         "[M-H2O3P]+",           "[M-HO3P+H]+",
    "[M-C2H4NO2]+",      "[M-C3H9N2]+",          "[M-C3H7NO]+",          "[M-C2H4NO2+H]+",       "[M-C3H6NO]+",
    "[M-4H2O+H]+",       "[M-C5H9]+",            "[M-C6H6N2O+CH3COO]+",  "[M-C2H5NO]+",          "[M-CH3N2O]+",
    "[M-C2H3O2]+",       "[M-C3H9N+H]+",         "[M-CH4N3]+",           "[M-C2H4NO]+",          "[M-C3H5O]+",
    "[2M-3H2O+H]+",      "[M-3H2O+H]+",          "[M-3H2O+2H]2+",        "[M-CH3S]+",            "[M-CH2O2]+",
    "[M-C2H5O]+",        "[M-CHO2]+",            "[M-HCOO]+",            "[M-CH2NO]+",           "[M-CH2NO+H]+",
    "[M-C2H4NO2+K]+",    "[2M-2H2O+H]+",         "[M-2H2O+H]+",          "[M-2H2O+2H]2+",        "[M-C3H6NO+K]+",
    "[M-CH5N]+",         "[M-CH3O]+",            "[M-C2H6N]+",           "[M-C2H6NO]+",          "[M-CH4N]+",
    "[M-C2H3O2+Cl]+",    "[M-HO3P+CH3COO]+",     "[M-C2H4NO+K]+",        "[M-H2O]+",             "[M-2H2O+NH4]+",
    "[M-H3N]+",          "[2M-H2O+H]+",          "[M-OH]+",              "[M-CH4]+",             "[M-NH2]+",
    "[M-HO+H]+",         "[M-CH3]+",             "[M-C3H7O+HCOO]+",      "[M-CH3+H]+",           "[M-CH3S+Cl]+",
    "[M-CHO2+K]+",       "[M-CH3N2O+CH3COO]+",   "[M]+",                 "[M]2+",                "[2M+H]+",
    "[M+H]+",            "[M+H]2+",              "[M+2H]+",              "[M+2H]2+",             "[M-CH3+NH4]+",
    "[M+3H]3+",          "[M-HO+Na]+",           "[M-H+Li]+",            "[M+Li]+",              "[M-CH3O+K]+",
    "[M-CH3+Na]+",       "[M-CH3O+HCOO]+",       "[M+CH3]+",             "[M+HO]+",              "[2M+NH4]+",
    "[3M+NH4]+",         "[M+NH4]+",             "[M-H2N+Cl]+",          "[M+H2O+H]+",           "[M+H+NH4]2+",
    "[M-HO+K]+",         "[M-H+Na]+",            "[2M+Na]+",             "[M+Na]+",              "[M+H+Na]+",
    "[M+H+Na]2+",        "[M+2H+Na]3+",          "[M-HO+HCOO]+",         "[M+CH3OH+H]+",         "[M-CHO2+Br]+",
    "[M+H2O+NH4]+",      "[2M+Ca-H]+",           "[3M+Ca-H]+",           "[M+Ca-H]+",            "[2M+K]+",
    "[3M+K]+",           "[M+K]+",               "[2M+Ca]2+",            "[3M+Ca]2+",            "[M+Ca]2+",
    "[M+H+K]2+",         "[M+H2O+Na]+",          "[2M+CH3CN+H]+",        "[M+CH3CN+H]+",         "[M-C2H3O]+",
    "[M+CH3CN+2H]2+",    "[2M-H+2Na]+",          "[M-H+2Na]+",           "[M+2Na]+",             "[M+2Na]2+",
    "[M+HCOO+H]+",       "[M+H+2Na]3+",          "[M+CH3OH+NH4]+",       "[M+CH3OH+Na]+",        "[2M+3H2O+2H]+",
    "[M+H2O+K]+",        "[M+CH3COOH+H]+",       "[M+C3H8O+H]+",         "[M-HO+Br]+",           "[M+HCOO+NH4]+",
    "[2M+CH3CN+Na]+",    "[M+CH3CN+Na]+",        "[M+HCOO+Na]+",         "[M+3Na]+",             "[M+3Na]3+",
    "[M+CH3OH+K]+",      "[M+CH3CN+NH4]+",       "[M+2K-H]+",            "[M+CH3COO+NH4]+",      "[M+C3H8O+NH4]+",
    "[M+C2H6OS+H]+",     "[M+CH3CN+K]+",         "[M+CH3COOH+Na]+",      "[M+C3H8O+Na]+",        "[M+2CH3CN+H]+",
    "[M+HCOO+K]+",       "[M+C3H8O+Na+H]+",      "[M+2CH3CN+2H]2+",      "[M+C2H6OS+NH4]+",      "[M+CH3COOH+K]+",
    "[M+C2H6OS+Na]+",    "[M+2CH3CN+2H]+",       "[M+3CH3CN+2H]2+",      "[M+3CH3CN+2H]+",
]


ION_NEG = [
    "[M-C6H3I2O]-",          "[M-C12H21O10]-",       "[M-C6H10O9P]-",         "[M-C5H13NO4P]-",       "[M-C6H8O6-H]-",
    "[M-C7H6N5O]-",          "[M-C7H5O5]-",          "[M-C6H11O5]-",          "[M-H3O6P2]-",          "[M-C7H5O4]-",
    "[M-C6H12O4]-",          "[M-C6H11O4]-",         "[M-C8H13O2]-",          "[M-C10H19]-",          "[M-C5H9O4]-",
    "[M-C5H8NO3]-",          "[M-I]-",               "[M-C6H7N2O]-",          "[M-C7H7O2]-",          "[M-C8H7O]-",
    "[M-C5H9O3]-",           "[M-C8H17]-",           "[M-CH4O4P]-",           "[M-C6H5O2]-",          "[M-C7H7O]-",
    "[M-C7H13]-",            "[M-H2O4P]-",           "[M-C7H11]-",            "[M-C3H6NO2]-",         "[M-H-CO2-2HF]-",
    "[M-H2O3S]-",            "[M-H2O3P]-",           "[M-C3H7O2]-",           "[M-C2H4NO2]-",         "[M-C3H6NO]-",
    "[M-C5H11]-",            "[M-C3H3O2]-",          "[M-C2H4O2]-",           "[M-C2H3O2]-",          "[M-CH2O2]-",
    "[M-C2H5O]-",            "[M-CHO2]-",            "[M-HCOO]-",             "[M-CH2NO]-",           "[M-CH3N2]-",
    "[M-C2H3O]-",            "[M-CH3O]-",            "[M-CHO]-",              "[M-Na]-",              "[M-H20-H]-",
    "[M-H-H2O]-",            "[M-H-NH3]-",           "[M-H2O]-",              "[M-OH]-",              "[M-CH4]-",
    "[M-NH2]-",              "[M-CH3]-",             "[M-3H]3-",              "[M-2H]-",              "[M-2H]2-",
    "[2M-H]-",               "[3M-H]-",              "[M-H]-",                "[M-H]2-",              "[M]-",
    "[M+HO]-",               "[2M-2H+Na]-",          "[M+Na-2H]-",            "[M+CH3OH-H]-",         "[2M+Cl]-",
    "[M+Cl]-",               "[M+H2O+H2O]-",         "[M+K-2H]-",             "[2M+HCOO-H]-",         "[M+HCOO-H]-",
    "[2M+HCOO]-",            "[M+HCOO]-",            "[M+CH3CN-H]-",          "[2M+CH3COO]-",         "[M+CH3COO]-",
    "[M+C3H8O-H]-",          "[M+CH3CN+Na-2H]-",     "[M+CH3OH+CH3OH]-",      "[M+C2H6OS-H]-",        "[2M+Br]-",
    "[M+Br]-",               "[M+CH3COONa-H]-",      "[M+HCOO+HCOO]-",        "[M+CH3COOH+CH3COOH]-", "[M+CH3CN+CH3CN]-",
    "[M+CF3COOH-H]-",        "[M+C3H8O+C3H8O]-",     "[M+C2H6OS+C2H6OS]-",
]


def parse_mf(formula: str) -> Dict[str, int]:
    """
    Parse a simple molecular formula (without parentheses) into a dict of element counts.

    Examples
    --------
    "C9H15"   -> {"C": 9, "H": 15}
    "C9H15M2" -> {"C": 9, "H": 15, "M": 2}
    "M"       -> {"M": 1}
    """
    pattern = re.compile(r'([A-Z][a-z]?)(\d*)')
    counts = Counter()

    for elem, num in pattern.findall(formula):
        n = int(num) if num else 1
        counts[elem] += n

    return dict(counts)


def calc_exact_mass(formula: str) -> float:
    """
    Calculate the exact monoisotopic mass of a formula string using EXACT_MASS.

    Parameters
    ----------
    formula : str
        Molecular formula string, e.g. "C9H15", "H2O", "C9H15M2".
        All element symbols in the formula must exist in EXACT_MASS.

    Returns
    -------
    float
        Exact mass (in unified atomic mass units, u).

    Notes
    -----
    - This implementation assumes simple formulas without parentheses,
      hydrates, or dot notation (e.g. "C6H12O6·H2O" is not supported).
    - If an element symbol is not found in EXACT_MASS,
      a ValueError is raised.
    """
    atoms = parse_mf(formula)
    total_mass = 0.0

    for elem, count in atoms.items():
        if elem not in EXACT_MASS:
            raise ValueError(
                f"Element {elem} is not found in EXACT_MASS. "
                f"Please add its exact mass before using it."
            )
        mass_elem = EXACT_MASS[elem]
        total_mass += mass_elem * count

    return total_mass


def calc_ion_mass(ion: str, M: float) -> float:
    """
    Calculate the exact monoisotopic mass of an ion expression.

    The ion expression is expected to be in square brackets with a
    charge at the end, for example:
        "[M+H]+"
        "[M-H]-"
        "[M+2H]2+"
        "[M-5H2O+H]+"
        "[M-2H]2-"
        "[2M-H]-"

    Rules (examples)
    ----------------
    Let M be a neutral molecule with known exact mass `M`, and
    H, H2O, etc. be fragments with known exact masses.

    - [M+H]+      : mass = M + H - e
    - [M-H]-      : mass = M - H + e
    - [M+2H]2+    : mass = M + 2*H - 2*e
    - [M-5H2O+H]+ : mass = M - 5*H2O + H - e
    - [M-2H]2-    : mass = M - 2*H + 2*e
    - [2M-H]-     : mass = 2*M - H + e

    Requirements
    ------------
    - EXACT_MASS must contain "e" as the electron mass (in u).
    - calc_exact_mass(formula: str) must be defined and able
      to handle formulas like "H", "H2O", etc.（不要求能处理前导系数，
      系数在本函数中拆解）
    - The symbol "M" is treated specially and uses the value `M`
      passed into this function (it is NOT taken from EXACT_MASS).

    Parameters
    ----------
    ion : str
        Ion expression, e.g. "[M+H]+", "[M-H]-", "[M+2H]2+", "[M-5H2O+H]+",
        "[M-2H]2-", "[2M-H]-", etc.
    M : float
        Exact mass of the neutral molecule M (in u).

    Returns
    -------
    float
        Exact mass of the ion (in unified atomic mass units, u).

    Raises
    ------
    ValueError
        If the format is invalid or charge cannot be parsed.
    """
    if "e" not in EXACT_MASS:
        raise ValueError("EXACT_MASS must contain the electron mass under key 'e'.")

    # 1) Extract the inner expression and the charge part
    match = re.match(r'^\[(.+)\]([+-]\d*|\d+[+-])$', ion.strip())
    if not match:
        raise ValueError(
            f"Invalid ion format: {ion!r}. "
            f"Expected something like '[M+H]+' or '[M+2H]2+'."
        )

    inner = match.group(1)      # e.g. "M+H", "M-5H2O+H", "M-2H", "2M-H"
    charge_str = match.group(2) # e.g. "+", "2+", "-2", "+3"

    # 2) Parse the charge value and sign
    # Allow formats: "+", "-", "+2", "-2", "2+", "3-"
    charge_match = re.match(r'^([+-])(\d*)$|^(\d+)([+-])$', charge_str)
    if not charge_match:
        raise ValueError(f"Cannot parse charge from {charge_str!r}.")

    if charge_match.group(1):  # form like "+", "-2", "+3"
        sign = 1 if charge_match.group(1) == '+' else -1
        num = charge_match.group(2)
        z = sign * (int(num) if num else 1)
    else:  # form like "2+", "3-", etc.
        num = int(charge_match.group(3))
        sign = 1 if charge_match.group(4) == '+' else -1
        z = sign * num

    # 3) Parse inner expression like "M+H-5H2O+H", "M-2H" or "2M-H"
    expr = inner.strip()
    if not expr.startswith(("+", "-")):
        expr = "+" + expr

    # Tokens are like "+M", "-H", "-5H2O", "+2H", "+2M", etc.
    tokens = re.findall(r'([+-][^+-]+)', expr)

    total_neutral_mass = 0.0

    for tok in tokens:
        sign_char = tok[0]
        part = tok[1:].strip()  # e.g. "M", "H", "5H2O", "2H", "2M"

        factor = 1 if sign_char == "+" else -1

        # 3.1 拆分前导系数（对所有片段都尝试，包括可能的 "2M", "5H2O"）
        coef = 1
        frag = part
        coef_match = re.match(r"^(\d+)([A-Za-z].*)$", part)
        if coef_match:
            coef = int(coef_match.group(1))
            frag = coef_match.group(2)

        # 3.2 计算该片段的质量
        if frag == "M":
            # 可能是 "M" 或 "2M" 等
            part_mass = coef * M
        else:
            # 其它片段（可能是 "H", "H2O" 等），交给 calc_exact_mass
            # 注意此时前导系数已拆掉，calc_exact_mass 只看到正常分子式
            part_mass = coef * calc_exact_mass(frag)

        total_neutral_mass += factor * part_mass

    # 4) Electron mass correction: ion_mass = neutral_mass - z * m_e
    # z > 0: remove z electrons -> subtract z * m_e
    # z < 0: add |z| electrons  -> add |z| * m_e
    m_e = EXACT_MASS["e"]
    ion_mass = total_neutral_mass - z * m_e

    return ion_mass


def is_valid_ion_type(ion: str, mf: str) -> bool:
    """
    Given a neutral molecular formula `mf` and an ion expression `ion`,
    check whether the ion is compositionally reasonable in terms of atom counts.

    Logic
    -----
    - Only the lost fragments (those with '-' sign inside the brackets) are checked.
      The total atoms in these fragments must not exceed the available atoms
      from the neutral molecule(s).
    - Supports leading integer coefficients such as 2H2O, 3CH3, 5C2H4, etc.
    - Supports 2M, 3M, ...: the available atoms are scaled by the M coefficient.
    - Electron loss/gain is ignored for atom balance; only element atom counts
      are considered.

    Parameters
    ----------
    ion : str
        Ion expression, e.g. "[M-CHO2]-", "[M-2H2O+H]+", "[2M-3H2O+H]+", etc.
    mf : str
        Neutral molecular formula of M, e.g. "C10H12O3N".

    Returns
    -------
    bool
        True if the ion is compositionally reasonable, False otherwise.
    """
    # 1) Parse the neutral molecule formula of M
    base_atoms = parse_mf(mf)  # e.g. {'C': x, 'H': y, ...}

    # 2) Parse ion: extract the inner expression within brackets
    m = re.match(r'^\[(.+)\]([+-]\d*|\d+[+-])$', ion.strip())
    if not m:
        # Invalid format is treated as unreasonable
        return False

    inner = m.group(1)   # e.g. "M-CHO2", "M-2H2O+H", "2M-3H2O+H"

    # 3) Normalize to start with '+' or '-', then split into tokens
    expr = inner.strip()
    if not expr.startswith(("+", "-")):
        expr = "+" + expr

    tokens = re.findall(r'([+-][^+-]+)', expr)

    # 4) Determine the coefficient of M (default 1)
    #    e.g. "M-..."  -> m_coeff = 1
    #         "2M-..." -> m_coeff = 2
    m_coeff = 0
    for tok in tokens:
        sign_char = tok[0]
        part = tok[1:].strip()  # "M", "2M", "H2O", "2H2O", ...

        if sign_char == "+":
            # Only '+' part contributes to total M coefficient
            coef = 1
            frag = part
            coef_match = re.match(r'^(\d+)([A-Za-z].*)$', part)
            if coef_match:
                coef = int(coef_match.group(1))
                frag = coef_match.group(2)

            if frag == "M":
                m_coeff += coef

    if m_coeff == 0:
        # If there is no M in the expression, we cannot compare to mf -> treat as invalid
        return False

    # 5) Accumulate all "lost fragments" (sign '-') atom counts
    lost_atoms = Counter()

    for tok in tokens:
        sign_char = tok[0]
        if sign_char != "-":
            continue

        part = tok[1:].strip()  # e.g. "CHO2", "2H2O", "3C2H4", "Na", ...

        # Try to parse leading integer coefficient, e.g. 2H2O
        coef = 1
        frag = part
        coef_match = re.match(r'^(\d+)([A-Za-z].*)$', part)
        if coef_match:
            coef = int(coef_match.group(1))
            frag = coef_match.group(2)

        # Try to parse frag as a normal molecular formula
        try:
            frag_atoms = parse_mf(frag)
        except Exception:
            # If the fragment cannot be parsed as a molecular formula, treat as invalid
            return False

        # Accumulate atoms in lost_atoms
        for elem, cnt in frag_atoms.items():
            lost_atoms[elem] += coef * cnt

    # 6) Compute the total available atoms, scaled by M coefficient
    #    For example, if the expression involves 2M, available atoms = base_atoms * 2
    total_available = {elem: cnt * m_coeff for elem, cnt in base_atoms.items()}

    # 7) Compare lost atoms vs. available atoms
    for elem, lost_cnt in lost_atoms.items():
        avail_cnt = total_available.get(elem, 0)
        if lost_cnt > avail_cnt:
            return False

    return True


def afford_ion_type(formula: str, ionmode: str) -> dict[str, float]:
    """
    For a given neutral molecular formula `mf` and ionization mode `ionmode`,
    generate all compositionally reasonable ion types and their exact masses.

    Parameters
    ----------
    mf : str
        Neutral molecular formula of M, e.g. "C10H12O3N".
    ionmode : str
        Ionization mode: "pos" for positive ions, "neg" for negative ions.

    Returns
    -------
    dict[str, float]
        A dictionary mapping ion type string to its exact monoisotopic mass.
        Example: { "[M+H]+": 123.0456, "[M+Na]+": 145.0278, ... }

    Notes
    -----
    - Ion types are taken from ION_POS or ION_NEG depending on ionmode.
    - Ion types are first filtered by `is_valid_ion_type(ion, mf)` so that the
      lost fragments do not require more atoms than available in M (or nM).
    - Neutral exact mass M is computed from `mf` using `calc_exact_mass`.
    """
    # 1) Choose the ion list according to ionmode
    if ionmode.lower() in ('positive', "pos", 'p', '+'):
        ion_list = ION_POS
    elif ionmode.lower() in ('negative', 'neg', 'n', '-'):
        ion_list = ION_NEG
    else:
        raise ValueError(f"Unsupported ionmode: {ionmode!r}. Use 'pos' or 'neg'.")

    # 2) Compute exact mass of neutral M
    M = calc_exact_mass(formula)

    # 3) Filter ion types by composition feasibility and compute ion masses
    result: dict[str, float] = {}

    for ion in ion_list:
        # Check if this ion type is compositionally reasonable for given mf
        if not is_valid_ion_type(ion, formula):
            continue

        # Compute ion mass using the given ion expression and neutral mass M
        mass = calc_ion_mass(ion, M)
        result[ion] = mass

    return result


def annotate_precursor_type(peakframe: PeakFrame,
                            mf: str,
                            ionmode:str, 
                            mz_on = 'PRECURSORMZ',
                            type_on = 'PRECURSORTYPE',                           
                            tol: float = 6e-6) -> PeakFrame:
    """
    Annotate precursor ion types based on measured PRECURSORMZ and a dict of
    theoretical ion masses, using relative mass error tolerance.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame that must contain a 'PRECURSORMZ' column (float).
    ions : dict[str, float]
        Mapping {ion_type: exact_mass}, e.g. {'[M-H]-': 179.0550, ...}.
    tol : float, optional
        Relative mass tolerance (dimensionless). Default 6e-6
        corresponds to 6 ppm (6e-6 = 6 / 1e6).

    Returns
    -------
    pd.DataFrame
        A copy of df with an extra column 'PRECURSORTYPE' containing the
        matched ion type (or None if no match within tolerance).
    """
    df = peakframe.copy()
    df[type_on] = None

    ions = afford_ion_type(mf, ionmode=ionmode)

    ion_types = np.array(list(ions.keys()))
    ion_masses = np.array(list(ions.values()), dtype=float)

    def match_one(mz: float) -> str | None:
        # relative error: |mz - mass| / mass
        rel_err = np.abs(ion_masses - mz) / ion_masses
        idx = np.argmin(rel_err)
        if rel_err[idx] <= tol:
            return ion_types[idx]
        return None

    df[type_on] = df[mz_on].apply(match_one)
    
    return df