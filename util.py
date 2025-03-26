def get_ion_label_superscript(i: str, c: int) -> str:
    return i + to_superscript(f"+{c}")


# Function to convert charge to a superscript representation
def to_superscript(s):
    superscript_map = {
        "1": "¹",
        "2": "²",
        "3": "³",
        "4": "⁴",
        "5": "⁵",
        "6": "⁶",
        "7": "⁷",
        "8": "⁸",
        "9": "⁹",
        "0": "⁰",
        "+": "⁺",
    }
    return "".join(superscript_map.get(c, "") for c in str(s))


def to_subscript(s):
    subscript_map = {
        "1": "₁",
        "2": "₂",
        "3": "₃",
        "4": "₄",
        "5": "₅",
        "6": "₆",
        "7": "₇",
        "8": "₈",
        "9": "₉",
        "0": "₀",
        "+": "₊",
        "-": "₋",
    }
    return "".join(subscript_map.get(c, "") for c in str(s))


def generate_fragmentation_latex(peptide, forward_indices, reverse_indices):
    """
    Generates a LaTeX string for visualizing peptide fragmentation sites.

    :param peptide: The peptide sequence.
    :param forward_indices: A list of indices for forward fragment ion sites.
    :param reverse_indices: A list of indices for reverse fragment ion sites.
    :return: A LaTeX string representing the peptide with fragmentation sites.
    """
    # Initialize an empty list to store each part of the LaTeX string
    latex_parts = []

    # Process each amino acid in the peptide
    for i, aa in enumerate(peptide):
        # Add the amino acid to the LaTeX parts
        latex_parts.append(aa)

        # Check if this position is a forward fragmentation site
        if i in forward_indices:
            latex_parts.append(r"_{\rfloor}")
        # Check if this position is a reverse fragmentation site
        # Note: Reverse indices count from the end of the peptide
        if (len(peptide) - 1 - i) in reverse_indices:
            latex_parts.append(r"^{\lceil}")

    # Join all parts into a single string and return
    return "".join(latex_parts)


def get_ion_label(i: str, c: int) -> str:
    return "+" * c + i