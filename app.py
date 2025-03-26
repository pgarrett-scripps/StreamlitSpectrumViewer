import tempfile
from typing import List

import pandas as pd
import streamlit_permalink as stp
import streamlit as st
import peptacular as pt

import matplotlib as mpl

from app_input import get_all_inputs, get_ion_label
import constants
from color_util import get_color_dict
from plot_util import (
    generate_annonated_spectra_plotly,
    coverage_string,
    generate_fragment_plot,
    generate_fragment_plot_ion_type,
)


@st.cache_data
def get_cached_fragments(annotation:pt.ProFormaAnnotation,
                  is_monoisotopic: bool,
                  fragment_types: list[str],
                  charges: list[int],
                  isotopes: list[int],
                  losses: list[(str, float)],
                  immonium_ions: bool):
    fragmenter = pt.Fragmenter(annotation, is_monoisotopic)

    fragments = pt.fragment(
        sequence=annotation,
        ion_types=fragment_types,
        charges=charges,
        monoisotopic=is_monoisotopic,
        isotopes=isotopes,
        losses=losses,
    )

    if immonium_ions:
        fragments.extend(
            fragmenter.fragment(
                ion_types=["i"],
                charges=[1],
                isotopes=isotopes,
                losses=losses,
            )
        )

    return fragments

st.set_page_config(page_title="Spectra Viewer", page_icon=":eyeglasses:", layout="wide")


with st.sidebar:


    c1, c2 = st.columns([3, 2], vertical_alignment="center")
    c1.title("Spectra Viewer :eyeglasses:")
    with c2:
        stateful = stp.toggle("Stateful", True, key="stateful")

    st.caption(
        'A tool to visualize and annotate msms spectra. Ensure to click the "Apply" button to update the '
        "visualization."
    )
    st.caption("Made with [Peptacular](https://pypi.org/project/peptacular/).")

    with stp.form("params"):

        btn = stp.form_submit_button("Apply", use_container_width=True, type='primary')

        if not stateful:
            st.query_params.clear()
            st.query_params["stateful"] = "false"

        params = get_all_inputs(stateful)


# check if charges are in bounds
if params.min_charge > params.max_charge:
    st.error(f"Minimum charge: {params.min_charge} cannot be greater than maximum charge: {params.max_charge}")
    st.stop()

# if more than 10 charges error
if params.max_charge - params.min_charge > constants.MAX_CHARGE_STATES:
    st.error(f"Cannot have more than {constants.MAX_CHARGE_STATES} charge states")
    st.stop()

# check if ppm error and tolerance are in bounds
if params.mass_tolerance_type == "ppm":
    if params.mass_tolerance > constants.MAX_PPM_MASS_TOLERANCE:
        st.error(f"Mass tolerance cannot be greater than {constants.MAX_PPM_MASS_TOLERANCE} for ppm")
        st.stop()

if params.mass_tolerance_type == "th":
    if params.mass_tolerance > constants.MAX_TH_MASS_TOLERANCE:
        st.error(f"Mass tolerance cannot be greater than {constants.MAX_TH_MASS_TOLERANCE} for th")
        st.stop()

# try to parse sequence
try:
    annotation = pt.parse(params.sequence)
except Exception as err:
    st.error(f"Error parsing peptide sequence: {err}")
    st.stop()

# try to get mas of annotation
try:
    _ = pt.mass(annotation, monoisotopic=True, ion_type="p", charge=0)
except Exception as err:
    st.error(f"Error calculating peptide mass: {err}")
    st.stop()

# see if there is any ambiguity in the sequence
if annotation.contains_sequence_ambiguity():
    st.error("Sequence cannot contain ambiguity!")
    st.stop()

color_dict = get_color_dict(params.min_charge, params.max_charge)

# Show Sequence Info
st.header(params.sequence)
c1, c2, c3, c4 = st.columns(4)

c1.metric("Mass", round(pt.mass(annotation), 4))
c2.metric("Length", len(params.unmodified_sequence))
c3.metric("Peaks", len(params.spectra))

cols = st.columns(len(params.charges))
for i in range(params.min_charge, params.max_charge + 1):
    col = cols[i - params.min_charge]
    col.metric(f"M/Z +{i}", round(pt.mz(params.sequence, charge=i), 4))

fragments = get_cached_fragments(annotation,
                          params.is_monoisotopic,
                          params.fragment_types,
                          params.charges,
                          params.isotopes,
                          params.losses,
                          params.immonium_ions)

c4.metric("Fragments", len(fragments))

frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

if params.spectra:

    mzs, ints = zip(*params.spectra)

    # Take top n peaks and bottom n peaks
    # top_n_spectra = sorted(zip(mzs, ints), key=lambda x: x[1], reverse=True)[:top_n]
    # bottom_n_spectra = sorted(zip(mzs, ints), key=lambda x: x[1])[:bottom_n]

    # Combine ensuring no duplicates - since it's a list of tuples, we can use a dictionary to remove duplicates
    # efficiently
    # spectra_dict = {mz: intensity for mz, intensity in top_n_spectra + bottom_n_spectra}
    spectra_dict = {mz: intensity for mz, intensity in zip(mzs, ints)}
    spectra = list(spectra_dict.items())
    mzs, ints = zip(*spectra)

    max_spectra_mz = max(mzs)

    # TODO: Add priority to fragment matches, a random isotope match should not be better than a non-isotope match
    fragment_matches = pt.get_fragment_matches(
        fragments,
        mzs,
        ints,
        params.mass_tolerance,
        params.mass_tolerance_type,
        "largest" if params.peak_assignment == "most intense" else "closest",
    )
    fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)

    if params.filter_missing_mono:
        fragment_matches = pt.filter_missing_mono_isotope(fragment_matches)

    if params.filter_interrupted_iso:
        fragment_matches = pt.filter_skipped_isotopes(fragment_matches)

    match_cov = pt.get_match_coverage(fragment_matches)

    if len(fragment_matches) == 0:
        st.warning(
            "No matches found, try increasing the mass tolerance or changing the ion types and charges"
        )
        st.stop()

    fragment_matches = {
        fm.mz: fm for fm in fragment_matches
    }  # keep the best fragment match for each mz

    match_data = []
    data = []
    for mz, i in zip(mzs, ints):
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

    match_df = pd.DataFrame(match_data)
    match_df["matched"] = True
    match_df["abs_error"] = match_df["error"].abs()
    match_df["abs_error_ppm"] = match_df["error_ppm"].abs()

    # spectra_df = spectra_df[~spectra_df['mz'].isin(match_df['mz'])]
    spectra_df = pd.concat([spectra_df, match_df])

    if params.hide_unassigned_peaks:
        spectra_df = spectra_df[spectra_df["matched"]]

    spectra_df["ion_color_type"] = spectra_df["ion_type"]
    spectra_df.loc[spectra_df["internal"], "ion_color_type"] = "i"

    def create_labels(row):
        if row["ion_type"] != "":
            charge_str = "+" * int(row["charge"])
            charge = int(row["charge"])
            ion_type_str = row["ion_type"]

            if row["internal"]:
                color_label = get_ion_label(row["ion_type"], int(row["charge"]))
                # ion_label = f"{get_ion_label(row['ion_type'], int(row['charge']))}{int(row['parent_number'])}i"
                ion_label = row["label"]
            else:
                color_label = get_ion_label(row["ion_type"], int(row["charge"]))
                # ion_label = f"{get_ion_label(row['ion_type'], int(row['charge']))}{int(row['parent_number'])}"
                ion_label = row["label"]

        else:
            color_label = "unassigned"
            ion_label = "unassigned"

        return ion_label, color_label

    # List comprehension to create ion and color labels
    labels = [create_labels(row) for _, row in spectra_df.iterrows()]
    spectra_df["ion_label"], spectra_df["color_label"] = zip(*labels)

    spectra_df.loc[spectra_df["ion_type"] == "I", "label"] = spectra_df.loc[
        spectra_df["ion_type"] == "I", "sequence"
    ].values

    # Assigning colors based on color labels
    spectra_df["color"] = [color_dict[label] for label in spectra_df["color_label"]]

    # change internal ions color to magenta
    # spectra_df.loc[spectra_df['internal'], 'color'] = 'magenta'

    cmap = mpl.colormaps.get_cmap("Blues")

    fig = generate_annonated_spectra_plotly(spectra_df, scale=params.y_axis_scale)
    st.plotly_chart(fig, use_container_width=True)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".svg") as tmpfile:
        # Save the figure to the temporary file
        fig.write_image(
            file=tmpfile.name, format="svg", width=1920, height=1080, scale=3.0
        )

        # Read the content of the temporary file
        tmpfile.seek(0)  # Go to the start of the file
        data = tmpfile.read()

    # Create a download button in Streamlit
    st.download_button(
        label="Download chart as SVG",
        data=data,
        file_name="spectra.svg",
        mime="image/svg+xml",
    )

    # Create a color map for the intensities
    st.caption("Sequence coverage")
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
                c = color_dict[get_ion_label(ion, charge)]
                s = coverage_string(cov_arr, params.unmodified_sequence, c)

                # center text
                ion_span = f'<span style="color:{c}">{ion}<sup>+{charge}</sup></span>'
                st.markdown(f"{ion_span} {s}", unsafe_allow_html=True)

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

            combined_data[get_ion_label(ion, charge)] = frags

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

    styled_dfs = []

    def color_by_ion_type(col):
        ion_type = col.name[-1]
        color = color_dict.get(
            ion_type, "grey"
        )  # get color or default to grey if not found
        return ["color: %s" % color] * len(col)

    for df in dfs:
        styled_df = df.style.apply(color_by_ion_type)

        # Set table styles with increased horizontal padding for more space between columns,
        # centered text, and no borders
        styles = [
            dict(
                selector="td",
                props=[
                    ("padding", "2px 2px"),
                    ("text-align", "center"),
                    ("border", "none"),
                ],
            ),
            dict(
                selector="th",
                props=[
                    ("padding", "2px 2px"),
                    ("text-align", "center"),
                    ("border", "none"),
                ],
            ),
        ]
        styled_df = styled_df.set_table_styles(styles)
        styled_dfs.append(styled_df)

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
                if ion in "abc":
                    ion_number = row + 1
                else:
                    ion_number = len(params.unmodified_sequence) - row
                label = col + str(ion_number)
                mz = data.loc[row, col]

                if mz <= params.min_mz or mz >= params.max_mz:
                    styled.loc[row, col] = (
                        f"background-color: #BEBEBE; color: black; text-align: center; font-weight: bold;"
                    )
                else:
                    if label in accepted_normal_ions:
                        styled.loc[row, col] = (
                            f"background-color: {color_dict[col]}; color: white; text-align: center; font-weight: bold;"
                        )
                    elif label in accepted_internal_ions:
                        styled.loc[row, col] = (
                            f"background-color: {color_dict[col]}; color: magenta; text-align: center; font-style: italic; font-weight: bold;"
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

    st.markdown("---")

    st.subheader("Fragment Ions")

    st.dataframe(combined_df, height=int(35.2 * (len_combined_df + 1)), hide_index=True)

    st.plotly_chart(
        generate_fragment_plot_ion_type(params.unmodified_sequence, spectra_df),
        use_container_width=True,
    )

    st.markdown("---")

    st.subheader("Stats")

    # calculate percentage if intensity accounted for
    total_intensity = spectra_df["intensity"].sum()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(label="Total Intensity", value=round(total_intensity, 1))
    c2.metric(label="Matched Intensity", value=round(match_df["intensity"].sum(), 1))
    c3.metric(
        label="Unmatched Intensity",
        value=round(spectra_df["intensity"].sum() - match_df["intensity"].sum(), 1),
    )
    c4.metric(
        label="Matched Intensity %",
        value=round(match_df["intensity"].sum() / total_intensity * 100, 2),
    )

    st.markdown("---")

    # error_fig = generate_error_histogram(spectra_df, mass_tolerance_type)
    # st.plotly_chart(error_fig, use_container_width=True)

    with st.expander("Fragments"):
        st.dataframe(frag_df)

    with st.expander("Peaks"):
        st.dataframe(spectra_df)


else:
    st.warning("No spectra....")
