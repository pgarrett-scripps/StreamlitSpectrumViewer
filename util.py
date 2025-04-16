import peptacular as pt
import pandas as pd

from app_input import SpectraInputs
from plot_util import coverage_string
import streamlit as st
from urllib.parse import quote_plus
import requests

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
        ".": "‧",
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

def get_ion_label_super(i: str, c: int) -> str:
    return f"<sup>+{c}</sup>{i}"

def get_fragment_matches(params: SpectraInputs, fragments: list[pt.Fragment]) -> list[pt.FragmentMatch]:
    mzs, ints = params.mz_int_values

    # Take top n peaks and bottom n peaks
    # top_n_spectra = sorted(zip(mzs, ints), key=lambda x: x[1], reverse=True)[:top_n]
    # bottom_n_spectra = sorted(zip(mzs, ints), key=lambda x: x[1])[:bottom_n]

    # Combine ensuring no duplicates - since it's a list of tuples, we can use a dictionary to remove duplicates
    # efficiently
    # spectra_dict = {mz: intensity for mz, intensity in top_n_spectra + bottom_n_spectra}
    spectra_dict = {mz: intensity for mz, intensity in zip(mzs, ints)}

    # TODO: Add priority to fragment matches, a random isotope match should not be better than a non-isotope match
    fragment_matches = pt.get_fragment_matches(
        fragments,
        mzs,
        ints,
        params.mass_tolerance,
        params.mass_tolerance_type,
        params.peak_assignment_type,
    )
    fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)

    if params.filter_missing_mono:
        fragment_matches = pt.filter_missing_mono_isotope(fragment_matches)

    if params.filter_interrupted_iso:
        fragment_matches = pt.filter_skipped_isotopes(fragment_matches)

    return fragment_matches


def get_match_cov(fragment_matches: list[pt.FragmentMatch]):
    return pt.get_match_coverage(fragment_matches)


def get_spectra_df(params: SpectraInputs, fragment_matches: list[pt.FragmentMatch]) -> pd.DataFrame:

    fragment_matches = {
        fm.mz: fm for fm in fragment_matches
    }  # keep the best fragment match for each mz

    match_data, data = [], []
    for mz, i in params.spectra:
        fm = fragment_matches.get(mz, None)
        if fm:
            match_data.append(fm.to_dict())
        else:
            fm = pt.FragmentMatch(fragment=None, mz=mz, intensity=i)
            d = fm.to_dict()
            
            # update any vlaues that are 0 to be None
            for key in d:
                if key not in ["mz", "intensity"]:
                    d[key] = None

            data.append(d)

    spectra_df = pd.DataFrame(data)
    spectra_df["matched"] = False
    spectra_df["abs_error"] = None
    spectra_df["abs_error_ppm"] = None

    if len(match_data) > 0:
        match_df = pd.DataFrame(match_data)
        match_df["matched"] = True
        match_df["abs_error"] = match_df["error"].abs()
        match_df["abs_error_ppm"] = match_df["error_ppm"].abs()
    else:
        match_df = pd.DataFrame()

    if params.hide_unassigned_peaks:
        spectra_df = spectra_df[spectra_df["matched"]]

    def create_label(row):
        # {charge}{ion_type}{number}{[isotope]}{(loss)}
        return f"{row['charge']}{row['ion_type']}{row['number']}" + (
            f"[{row['isotope']}]" if row["isotope"] != 0 else ""
        ) + (f"({row['loss']})" if row["loss"] != 0 else "")


    match_df['ion_label'] = match_df.apply(create_label, axis=1)
    spectra_df['ion_label'] = ''

    def create_ion_group_label(row):
        # {charge}{ion_type}
        return f"{row['charge']}{row['ion_type']}"
    
    match_df['ion_group_label'] = match_df.apply(create_ion_group_label, axis=1)
    spectra_df['ion_group_label'] = 'unassigned'

    grouped_df = pd.concat([spectra_df, match_df])

    # Assigning colors based on color labels
    grouped_df["color"] = [params.get_color(i, c) for i, c in zip(grouped_df["ion_type"], grouped_df["charge"])]

    return grouped_df


def get_spectra_dfold(params: SpectraInputs, fragment_matches: list[pt.FragmentMatch]) -> pd.DataFrame:

    fragment_matches = {
        fm.mz: fm for fm in fragment_matches
    }  # keep the best fragment match for each mz

    match_data, data = [], []
    for mz, i in params.spectra:
        fm = fragment_matches.get(mz, None)
        if fm:
            match_data.append(fm.to_dict())
        else:
            fm = pt.FragmentMatch(fragment=None, mz=mz, intensity=i)
            data.append(fm.to_dict())

    spectra_df = pd.DataFrame(data)
    spectra_df["matched"] = False
    spectra_df["abs_error"] = spectra_df["error"].abs()
    spectra_df["abs_error_ppm"] = spectra_df["error_ppm"].abs()

    if len(match_data) > 0:
        match_df = pd.DataFrame(match_data)
        match_df["matched"] = True
        match_df["abs_error"] = match_df["error"].abs()
        match_df["abs_error_ppm"] = match_df["error_ppm"].abs()
    else:
        match_df = pd.DataFrame()

    spectra_df = pd.concat([spectra_df, match_df])

    if params.hide_unassigned_peaks:
        spectra_df = spectra_df[spectra_df["matched"]]


    spectra_df["ion_color_type"] = spectra_df["ion_type"]
    spectra_df.loc[spectra_df["internal"], "ion_color_type"] = "i"
    

    def create_labels(row):
        if row["ion_type"] != "":
            color_label = get_ion_label(row["ion_type"], int(row["charge"]))
            ion_label = get_ion_label(row["ion_type"], int(row["charge"]))
        else:
            color_label = "unassigned"
            ion_label = "unassigned"

        return ion_label, color_label

    def create_labels_super(row):
        if row["ion_type"] != "":
            charge_str = "+" + str(row["charge"])
            charge = int(row["charge"])
            ion_type_str = row["ion_type"]
            ion_number = row["number"]
            isotope_number = row["isotope"]
            loss = row["loss"]

            # color_label = get_ion_label(row["ion_type"], int(row["charge"]))
            # ion_label = f"{get_ion_label(row['ion_type'], int(row['charge']))}{int(row['parent_number'])}"
            color_label = f"<sup>{charge_str}</sup>{ion_type_str}"
            ion_label = f"<sup>{charge_str}</sup>{ion_type_str}<sub>{ion_number}</sub>"

            if isotope_number != 0:
                ion_label += f"<sub>[{isotope_number}]</sub>"

            if loss != 0:
                ion_label += f"<sub>({str(round(row['loss'], 2))})</sub>"

            
        else:
            color_label = "unassigned"
            ion_label = "unassigned"

        return ion_label, color_label


    def create_labels_formated(row):
        # load pip install mathjax
        
        if row["ion_type"] != "":
            charge = int(row["charge"])
            ion_type = row["ion_type"]
            
            # Build LaTeX string components
            charge_latex = f"^{{{'+' + str(charge)}}}" if charge > 0 else ""
            
            # Check if isotope exists and has a valid value
            isotope_latex = ""
            if "isotope" in row and pd.notna(row["isotope"]) and row["isotope"] > 0:
                isotope_latex = f"^{{{'*'+str(row['isotope'])}}}"
            
            # Check if loss exists and has a valid value
            loss_latex = ""
            if "loss" in row and pd.notna(row["loss"]) and row["loss"] != 0:
                loss_latex = f"_{{({str(round(row['loss'],2))})}}"
                
            ion_number_latex = f"_{{{'i' if row['internal'] else ''}{row['number']}}}"
            
            # Create formatted ion label in LaTeX
            ion_label = f"${isotope_latex}{loss_latex}{ion_type}{charge_latex}{ion_number_latex}$"
            color_label = f"${ion_type}{charge_latex}$"
        else:
            color_label = "unassigned"
            ion_label = "unassigned"

        return ion_label, color_label

    # List comprehension to create ion and color labels
    labels = [create_labels(row) for _, row in spectra_df.iterrows()]
    spectra_df["ion_label"], spectra_df["color_label"] = zip(*labels)

    formatted_labels = [
        create_labels_formated(row) for _, row in spectra_df.iterrows()
    ]
    spectra_df["ion_label_formated"], spectra_df["color_label_formated"] = zip(
        *formatted_labels
    )

    super_labels = [
        create_labels_super(row) for _, row in spectra_df.iterrows()
    ]
    spectra_df["ion_label_super"], spectra_df["color_label_super"] = zip(
        *super_labels
    )
    spectra_df["label"] = spectra_df["ion_label_super"]


    spectra_df.loc[spectra_df["ion_type"] == "I", "label"] = spectra_df.loc[
        spectra_df["ion_type"] == "I", "sequence"
    ].values

    # Assigning colors based on color labels
    spectra_df["color"] = [params.color_dict[label] for label in spectra_df["color_label"]]
    spectra_df["color_label"] = spectra_df["color_label_super"]

    return spectra_df


def display_coverage_markdown(params: SpectraInputs, spectra_df: pd.DataFrame):
    for ion in params.fragment_types:
        for charge in params.charges:
            cov_arr = [0] * len(params.unmodified_sequence)
            tmp_df = spectra_df[
                (spectra_df["ion_type"] == ion) & (spectra_df["charge"] == charge)
                ]

            if ion in "abc":
                for num in tmp_df["end"].unique():
                    cov_arr[num - 1] = 1
            elif ion in "xyz":
                for num in tmp_df["start"].unique():
                    cov_arr[num] = 1
            else:
                continue

            if len(cov_arr) > 0:
                c = params.color_dict[get_ion_label(ion, charge)]
                s = coverage_string(cov_arr, params.unmodified_sequence, c)

                # center text
                ion_span = f'<span style="color:{c}">{ion}<sup>+{charge}</sup></span>'
                st.markdown(f"{ion_span} {s}", unsafe_allow_html=True)


def get_fragment_match_table(params: SpectraInputs, spectra_df: pd.DataFrame, frag_df: pd.DataFrame) -> pd.DataFrame:
    dfs = []
    # combined_data = {'AA': list(unmodified_sequence)}
    combined_data = {"AA": pt.split(params.sequence)}
    for ion in params.fragment_types:
        for charge in params.charges:
            data = {"AA": pt.split(params.sequence)}
            ion_df = frag_df.copy()
            ion_df = ion_df[
                (ion_df["ion_type"] == ion)
                & (ion_df["charge"] == charge)
                & (ion_df["internal"] == False)
                & (ion_df["isotope"] == 0)
                & (ion_df["loss"] == 0)
                ]
            ion_df.sort_values(
                by=["start"] if ion in "xyz" else ["end"],
                inplace=True,
                ascending=False if ion in "xyz" else True,
            )

            # keep only a single number
            ion_df.drop_duplicates(
                subset=["start"] if ion in "xyz" else ["end"], inplace=True
            )

            frags = ion_df["mz"].tolist()

            if ion in "xyz":
                frags = frags[::-1]

            data[ion] = frags

            combined_data[f"{charge}{ion}"] = frags

            # Displaying the table
            df = pd.DataFrame(data)
            df["# (abc)"] = list(range(1, len(df) + 1))
            df["# (xyz)"] = list(range(1, len(df) + 1))[::-1]

            # reorder columns so that # is first # +1 is last and AA is in the middle
            df = df[
                ["AA"]
                + ["# (abc)"]
                + [col for col in df.columns if col not in ["AA", "# (abc)", "# (xyz)"]]
                + ["# (xyz)"]
                ]
            dfs.append(df)


    combined_df = pd.DataFrame(combined_data)
    # sort columns based on alphabetical order
    combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

    def highlight_cells(data):
        # Initialize empty DataFrame with same index and columns as original
        styled = pd.DataFrame("", index=data.index, columns=data.columns)

        # Iterate over cells and update `styled` based on cell position
        for row in data.index:
            for col in data.columns:
                if col == "AA" or col == "# (abc)" or col == "# (xyz)":
                    styled.loc[row, col] = (
                        f"background-color: gainsboro; color: black; text-align: center; font-weight: bold;"
                    )
                    continue

                ion = col[-1]
                charge = int(col[:-1])
                if ion in "abc":
                    ion_number = row + 1
                else:
                    ion_number = len(params.unmodified_sequence) - row
                ion_key = col + str(ion_number)
                mz = data.loc[row, col]

                if mz <= params.min_mz or mz >= params.max_mz:
                    styled.loc[row, col] = (
                        f"background-color: #BEBEBE; color: black; text-align: center; font-weight: bold;"
                    )
                else:
                    if ion_key in accepted_normal_ions:
                        styled.loc[row, col] = (
                            f"background-color: {params.get_color(ion, charge)}; color: white; text-align: center; font-weight: bold;"
                        )
                    elif ion_key in accepted_internal_ions:
                        styled.loc[row, col] = (
                            f"background-color: {params.get_color(ion, charge)}; color: magenta; text-align: center; font-style: italic; font-weight: bold;"
                        )
                    else:
                        styled.loc[row, col] = (
                            f"background-color: white; color: black; text-align: center;"
                        )
        return styled

    matched_ions = spectra_df[spectra_df["ion_type"] != ""]
    accepted_normal_ions = matched_ions[matched_ions["internal"] == False][
        "ion_label"
    ].tolist()
    accepted_internal_ions = matched_ions[matched_ions["internal"] == True][
        "ion_label"
    ].tolist()
    accepted_internal_ions = [ion[:-1] for ion in accepted_internal_ions]

    combined_df["# (abc)"] = list(range(1, len(params.unmodified_sequence) + 1))
    combined_df["# (xyz)"] = list(range(1, len(params.unmodified_sequence) + 1))[::-1]

    # reorder columns so that # is first # +1 is last and AA is in the middle
    combined_cols = combined_df.columns.tolist()
    combined_cols.remove("# (abc)")
    combined_cols.remove("# (xyz)")
    combined_cols.remove("AA")
    forward_cols = [
        col for col in combined_cols if "a" in col or "b" in col or "c" in col
    ]
    reverse_cols = [
        col for col in combined_cols if "x" in col or "y" in col or "z" in col
    ]

    # sort
    forward_cols.sort()
    reverse_cols.sort(reverse=True)

    new_cols = ["# (abc)"] + forward_cols + ["AA"] + reverse_cols + ["# (xyz)"]
    combined_df = combined_df[new_cols]
    len_combined_df = len(combined_df)

    combined_df = combined_df.style.format(precision=4).apply(
        highlight_cells, axis=None
    )

    return combined_df


def get_query_params_url(params_dict):
    """
    Create url params from alist of parameters and a dictionary with values.

    Args:
        params_list (str) :
            A list of parameters to get the value of from `params_dict`
        parmas_dict (dict) :
            A dict with values for the `parmas_list .
        **kwargs :
            Extra keyword args to add to the url
    """
    return "?" + "&".join(
        [
            f"{key}={quote_plus(str(value))}"
            for key, values in params_dict.items()
            for value in listify(values)
        ]
    )


def listify(o=None):
    if o is None:
        res = []
    elif isinstance(o, list):
        res = o
    elif isinstance(o, str):
        res = [o]
    else:
        res = [o]
    return res


def shorten_url(url: str) -> str:
    """Shorten a URL using TinyURL."""
    api_url = f"http://tinyurl.com/api-create.php?url={url}"

    try:
        response = requests.get(api_url)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        return f"Error: {e}"