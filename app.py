import tempfile

import pandas as pd
import streamlit_permalink as stp
import streamlit as st
import peptacular as pt

import matplotlib as mpl

from app_input import get_all_inputs
import constants
from color_util import get_color_dict
from plot_util import (
    generate_annonated_spectra_plotly,
    generate_fragment_plot_ion_type, generate_error_histogram,
)
from util import get_fragment_matches, get_match_cov, get_spectra_df, display_coverage_markdown, \
    get_fragment_match_table, get_query_params_url, shorten_url
from streamlit_js_eval import get_page_location


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

    st.markdown(f"""
            <div style='text-align: center; padding: 15px; top-margin: 0px'>
                <h3 style='margin: 0; font-size: 1.5em; color: #333;'>Spec-Viewer 👓</h3>
                <p style='font-size: 1.1em; line-height: 1.6; color: #555;'>
                    Powered by 
                    <a href="https://github.com/pgarrett-scripps/peptacular" target="_blank" style='color: #007BFF; text-decoration: none;'>
                        <strong>Peptacular</strong>
                    </a>. 
                    See the 
                    <a href="https://peptacular.readthedocs.io/en/latest/modules/getting_started.html#proforma-notation" 
                    target="_blank" style='color: #007BFF; text-decoration: none;'>
                        Proforma Notation Docs
                    </a> for supported peptide syntax. To report any issues or suggest improvements, please visit the 
                    <a href="https://github.com/pgarrett-scripps/StreamlitSpectrumViewer" 
                    target="_blank" style='color: #007BFF; text-decoration: none;'>
                        PepFrag Github Repo.
                    </a>
                </p>
            </div>
        """, unsafe_allow_html=True)

    stateful = stp.toggle("Stateful", True, key="stateful")

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

fragments = get_cached_fragments(annotation,
                          params.is_monoisotopic,
                          params.fragment_types,
                          params.charges,
                          params.isotopes,
                          params.losses,
                          params.immonium_ions)

frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

if not params.spectra:
    st.warning("No spectra....")
    st.stop()

fragment_matches = get_fragment_matches(params, fragments)

if not fragment_matches:
    st.warning(
        "No matches found, try increasing the mass tolerance or changing the ion types and charges"
    )

color_dict = get_color_dict(params.min_charge, params.max_charge)
match_cov = get_match_cov(fragment_matches)
spectra_df = get_spectra_df(params, fragment_matches)
match_df = spectra_df[spectra_df["matched"]]

cmap = mpl.colormaps.get_cmap("Blues")
spectra_fig = generate_annonated_spectra_plotly(spectra_df, scale=params.y_axis_scale,
                                                error_scale=params.mass_tolerance_type,
                                                line_width=params.line_width, 
                                                text_size=params.text_size,
                                                marker_size=params.marker_size)
combined_df = get_fragment_match_table(params, spectra_df, frag_df)

top_window, bottom_window = st.container(), st.container()

with bottom_window:
    page_loc = get_page_location()

with top_window:
    title_c, _, button_c = st.columns([2, 1, 1])
    title_c.header("SpecView Results")

    if params.stateful:
        st.caption(
            '''**This pages URL automatically updates with your input, and can be shared with others. 
           You can optionally use the Generate TinyURL button to create a shortened URL.**''',
            unsafe_allow_html=True,
        )

        if page_loc and 'origin' in page_loc:
            url_origin = page_loc['origin']

            if button_c.button("Generate TinyURL", key="generate_tinyurl", type="primary"):
                url_params = {k: st.query_params.get_all(k) for k in st.query_params.keys()}
                page_url = f"{url_origin}{get_query_params_url(url_params)}"
                short_url = shorten_url(page_url)


                @st.dialog(title="Share your results")
                def url_dialog(url):
                    st.write(f"Shortened URL: {url}")


                url_dialog(short_url)

    st.divider()

    # Show Sequence Info
    st.subheader(params.sequence)
    c1, c2, c3, c4 = st.columns(4)

    c1.metric("Mass", round(pt.mass(annotation), 4))
    c2.metric("Peaks", len(params.spectra))
    c3.metric("Fragments", len(fragments))
    total_intensity = spectra_df["intensity"].sum()
    #c1.metric(label="Total Intensity", value=round(total_intensity, 1))
    #c2.metric(label="Matched Intensity", value=round(match_df["intensity"].sum(), 1))
    #c3.metric(label="Unmatched Intensity", value=round(spectra_df["intensity"].sum() - match_df["intensity"].sum(), 1))
    c4.metric(
        label="Matched Intensity %",
        value=round(match_df["intensity"].sum() / total_intensity * 100, 2),
    )

    spectra_tab, coverage_tab, data_tab = st.tabs(["Spectra", "Coverage", "Data"])

    with spectra_tab:

        @st.fragment
        def run_spectra(spectra_fig, params):
            # slider to zoom from min to max mz
            min_mz, max_mz = params.min_spectra_mz, params.max_spectra_mz
            mz_range = st.slider("Zoom M/Z Range", min_value=min_mz, max_value=max_mz, value=(min_mz, max_mz),
                                 key="plot_mz_range")



            # update the spectra_df with the new mz range
            # update fid to be from 100-200 mz
            spectra_fig = spectra_fig.update_layout(
                xaxis=dict(range=[mz_range[0], mz_range[1]]),
                xaxis2=dict(range=[mz_range[0], mz_range[1]])  # If you have multiple x-axes
            )

            st.plotly_chart(spectra_fig, use_container_width=True)

        run_spectra(spectra_fig, params)

        with tempfile.NamedTemporaryFile(delete=False, suffix=".svg") as tmpfile:
            # Save the figure to the temporary file
            spectra_fig.write_image(
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
            use_container_width=True,
            on_click="ignore",
        )


    with coverage_tab:

        st.subheader("Sequence Coverage", divider=True)
        display_coverage_markdown(params, spectra_df)

        st.subheader("Fragment Matches", divider=True)

        df_html = combined_df.to_html()

        for col in combined_df.columns:
            if col.startswith("+"):
                df_html = df_html.replace(f'{col}</th>', f'{col.replace("+", "")}<sup>{col.count("+")}+</sup></th>')

        st.html(df_html)

        st.subheader("Fragment Locations", divider=True)
        st.plotly_chart(
            generate_fragment_plot_ion_type(params.unmodified_sequence, spectra_df),
            use_container_width=True,
        )

        #error_fig = generate_error_histogram(spectra_df, params.mass_tolerance_type)
        #st.plotly_chart(error_fig, use_container_width=True)

    with data_tab:

        st.subheader("Fragment Data", divider=True)
        st.dataframe(frag_df, hide_index=True)

        st.download_button(
            label="Download Data",
            data=frag_df.to_csv(index=False).encode("utf-8"),
            file_name=f"{annotation.serialize()}_fragment_data.csv",
            use_container_width=True,
            type="secondary",
            on_click="ignore",
            key="download_frag_data",
        )

        st.subheader("Spectra Data", divider=True)
        st.dataframe(spectra_df, hide_index=True)

        st.download_button(
            label="Download Data",
            data=spectra_df.to_csv(index=False).encode("utf-8"),
            file_name=f"{annotation.serialize()}_spectra_data.csv",
            use_container_width=True,
            type="secondary",
            on_click="ignore",
            key="download_spectra_data",
        )


    st.divider()

    st.markdown(f"""
        <div style='display: flex; justify-content: space-between; align-items: center; padding: 15px 0; border-top: 0px solid #ddd;'>
            <div style='text-align: left; font-size: 1.1em; color: #555;'>
                <a href="https://github.com/pgarrett-scripps/pep-frag" target="_blank" 
                   style='text-decoration: none; color: #007BFF; font-weight: bold;'>
                    Spec-Viewer
                </a>
                <a href="https://doi.org/10.5281/zenodo.15092689" target="_blank" style="margin-left: 12px;">
                    <img src="https://zenodo.org/badge/728447115.svg" alt="DOI" 
                         style="vertical-align: middle; height: 20px;">
                </a>
            </div>
            <div style='text-align: right; font-size: 1.1em; color: #555;'>
                <a href="https://github.com/pgarrett-scripps/peptacular" target="_blank" 
                   style='text-decoration: none; color: #007BFF; font-weight: bold;'>
                    Peptacular
                </a>
                <a href="https://doi.org/10.5281/zenodo.15054278" target="_blank" style="margin-left: 12px;">
                    <img src="https://zenodo.org/badge/591504879.svg" alt="DOI" 
                         style="vertical-align: middle; height: 20px;">
                </a>
            </div>
        </div>
    """, unsafe_allow_html=True)