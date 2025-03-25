# TODO: Add loss analysis plot

import tempfile
from typing import List

import pandas as pd
import streamlit_permalink as st
import peptacular as pt

import matplotlib as mpl
import ms_deisotope

import constants
from color_util import get_color_dict
from plot_util import generate_annonated_spectra_plotly, coverage_string, generate_fragment_plot, \
    generate_fragment_plot_ion_type
from query_params import parse_query_params, InvalidQueryParam, generate_app_url, QueryParams

st.set_page_config(page_title="Spectra Viewer", page_icon=":eyeglasses:", layout="wide")

try:
    qp = parse_query_params(st.query_params)
except InvalidQueryParam as e:
    st.error(str(e))
    st.stop()

# Initialize session state for fragment types if it doesn't exist (page refresh)
if 'fragment_types' not in st.session_state:
    st.session_state.fragment_types = set()
    st.session_state.fragment_types.update(qp.fragment_types)

if 'internal_fragment_types' not in st.session_state:
    st.session_state.internal_fragment_types = set()
    st.session_state.internal_fragment_types.update(qp.internal_fragment_types)


def update_fragment_types(label, is_checked):
    if is_checked:
        st.session_state.fragment_types.add(label)
    elif label in st.session_state.fragment_types:
        st.session_state.fragment_types.remove(label)


def update_internal_fragment_types(label, is_checked):
    if is_checked:
        st.session_state.internal_fragment_types.add(label)
    elif label in st.session_state.internal_fragment_types:
        st.session_state.internal_fragment_types.remove(label)


def get_ion_label_superscript(i: str, c: int) -> str:
    return i + to_superscript(f'+{c}')


# Function to convert charge to a superscript representation
def to_superscript(s):
    superscript_map = {
        "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵",
        "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "0": "⁰",
        '+': '⁺'
    }
    return ''.join(superscript_map.get(c, '') for c in str(s))


def to_subscript(s):
    subscript_map = {
        "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅",
        "6": "₆", "7": "₇", "8": "₈", "9": "₉", "0": "₀",
        '+': '₊', '-': '₋'
    }
    return ''.join(subscript_map.get(c, '') for c in str(s))


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
    return ''.join(latex_parts)


# Add query option for selected mass_tolerance to st.session_state. This fixes the issue when switching between ppm
# and th mass tolerance types, and causes the mass tolerance to be too large / too small
if qp.mass_tolerance_type == 'ppm' and 'ppm_mass_error' not in st.session_state:
    st.session_state['ppm_mass_error'] = qp.mass_tolerance
if qp.mass_tolerance_type == 'ppm' and 'th_mass_error' not in st.session_state:
    st.session_state['th_mass_error'] = 0.5
if qp.mass_tolerance_type == 'th' and 'th_mass_error' not in st.session_state:
    st.session_state['th_mass_error'] = qp.mass_tolerance
if qp.mass_tolerance_type == 'th' and 'ppm_mass_error' not in st.session_state:
    st.session_state['ppm_mass_error'] = 30.0


def get_ion_label(i: str, c: int) -> str:
    return '+' * c + i


with st.sidebar:
    st.title('Spectra Viewer :eyeglasses:')

    st.caption('A tool to visualize and annotate msms spectra. Ensure to click the "Apply" button to update the '
                'visualization.')

    st.caption('Made with [Peptacular](https://pypi.org/project/peptacular/).')

    sequence = st.text_input(label='Sequence (Proforma2.0 Notation)',
                             value=qp.sequence,
                             help=constants.SEQUENCE_HELP).replace(' ', '')
    unmodified_sequence = pt.strip_mods(sequence)

    c1, c2 = st.columns(2)
    mass_tolerance_type = c1.selectbox(
        label='Mass Tolerance Type',
        options=['ppm', 'th'],
        index=0 if qp.mass_tolerance_type == 'ppm' else 1,
        help=constants.MASS_TOLERANCE_TYPE_HELP
    )

    mass_tolerance = c2.number_input(
        label='Mass Tolerance',
        value=st.session_state['ppm_mass_error'] if mass_tolerance_type == 'ppm' else st.session_state['th_mass_error'],
        max_value=1000.0 if mass_tolerance_type == 'ppm' else 1.0,
        min_value=0.0,
        help=constants.MASS_TOLERANCE_HELP)

    # Update session state for mass tolerances
    if mass_tolerance_type == 'ppm':
        st.session_state['ppm_mass_error'] = mass_tolerance
    else:
        st.session_state['th_mass_error'] = mass_tolerance

    spectra = st.text_area(
        label='Spectra (mz Intensity)',
        value='\n'.join(f'{s[0]} {s[1]}' for s in qp.spectra),
        help=constants.SPECTRA_HELP,
        height=200,
    )

    st.subheader('Fragmentation Options')

    c1, c2 = st.columns(2)
    min_charge = c1.number_input(label='Min Charge',
                                 value=qp.min_charge,
                                 min_value=constants.MIN_CHARGE,
                                 max_value=constants.MAX_CHARGE)

    max_charge = c2.number_input(label='Max Charge',
                                 value=qp.max_charge,
                                 min_value=constants.MIN_CHARGE,
                                 max_value=constants.MAX_CHARGE)

    color_dict = get_color_dict(min_charge, max_charge)

    c1, c2 = st.columns(2)
    deselected_all = c1.button(label='Deselect All Ions', use_container_width=True)
    select_all = c2.button(label='Select All Ions', use_container_width=True)

    with st.form(key='query_form'):

        c1, c2 = st.columns(2)

        c1.subheader('Options Form')
        c2.form_submit_button(label='Apply', use_container_width=True, type='primary')

        st.write('Terminal Ions')
        fragment_types = set()
        for i in constants.IONS:
            cols = st.columns(max_charge - min_charge + 1)
            for c, col in enumerate(cols, min_charge):
                label = get_ion_label(i, c)
                val = label in st.session_state.fragment_types
                if deselected_all:
                    val = False
                if select_all:
                    val = True
                ion_selected = col.checkbox(label=get_ion_label_superscript(i, c),
                                            value=val,
                                            key=label)
                if ion_selected:
                    fragment_types.add(label)
                update_fragment_types(label, ion_selected)

        immonium_ions = st.checkbox(
            label='Immonium Ions',
            value=qp.immonium_ions,
            help=constants.IMMONIUM_IONS_HELP)

        internal_fragment_types = set()
        with st.expander('Internal Ions'):
            st.write('Internal Ions')
            for i in constants.INTERNAL_IONS:
                cols = st.columns(max_charge - min_charge + 1)
                for c, col in enumerate(cols, min_charge):
                    label = get_ion_label(i, c)
                    val = label in st.session_state.internal_fragment_types
                    if deselected_all:
                        val = False
                    ion_selected = col.checkbox(label=get_ion_label_superscript(i, c),
                                                value=val,
                                                key=label)
                    if ion_selected:
                        internal_fragment_types.add(label)

                    update_internal_fragment_types(label, ion_selected)

        internal_fragments = len(internal_fragment_types) > 0

        st.caption('Neutral Losses')
        cols = st.columns(len(constants.NEUTRAL_LOSSES))
        losses = []
        neutral_losses = []
        for i, (nl, mass) in enumerate(constants.NEUTRAL_LOSSES.items()):
            l = cols[i].checkbox(
                label=nl,
                value=nl in qp.neutral_losses,
                key=nl
            )
            if l:
                losses.append(mass)
                neutral_losses.append(nl)

        custom_losses = st.text_input(
            label='Custom Losses',
            value=';'.join([f'{nl[0]}:{nl[1]}' for nl in qp.custom_losses]),
            help=constants.CUSTOM_LOSSES_HELP
        )

        if custom_losses:
            custom_losses = [(nl.split(':')[0], float(nl.split(':')[1])) for nl in custom_losses.split(';')]
            losses.extend(custom_losses)

        c1, c2 = st.columns(2)
        mass_type = c1.selectbox(
            label='Mass Type',
            options=['monoisotopic', 'average'],
            index=0 if qp.mass_type == 'monoisotopic' else 1,
            help=constants.MASS_TYPE_HELP
        )

        peak_assignment = c2.selectbox(
            label='Peak Assignment',
            options=['largest', 'closest'],
            index=0 if qp.peak_assignment == 'largest' else 1,
            help=constants.PEAK_ASSIGNMENT_HELP)

        with st.expander('Isotopes'):

            isotopes = st.number_input(
                label='Isotopes',
                value=qp.isotopes,
                min_value=0,
                max_value=5,
                help=constants.ISOTOPES_HELP
            )

            c1, c2 = st.columns(2)
            filter_missing_mono = c1.checkbox(
                label='Filter missing mono peaks',
                value=qp.filter_missing_mono,
                help=constants.FILTER_MISSING_MONO_HELP
            )
            filter_interrupted_iso = c2.checkbox(
                label='Filter interrupted isotopes',
                value=qp.filter_interrupted_iso,
                help=constants.FILTER_INTERRUPTED_ISO_HELP
            )

        with st.expander('Plot Options'):

            c1, c2 = st.columns(2)

            y_axis_scale = c1.radio(
                label='Y Axis Scale',
                options=['linear', 'log'],
                horizontal=True,
                index=0 if qp.y_axis_scale == 'linear' else 1,
                help=constants.Y_AXIS_SCALE_HELP
            )

            hide_unassigned_peaks = c2.checkbox(
                label='Hide Unassigned Peaks',
                value=qp.hide_unassigned_peaks,
                help=constants.HIDE_UNASSIGNED_PEAKS_HELP
            )

        with st.expander('Peak Picker'):

            peak_picker = st.checkbox(
                label='Peak Picker',
                value=qp.peak_picker,
                help=constants.PEAK_PICKER_HELP)

            c1, c2 = st.columns(2)
            peak_picker_min_intensity = c1.number_input(
                label='Peak Picker Min Intensity',
                value=qp.peak_picker_min_intensity,
                min_value=0.0,
                max_value=1e9,
                help=constants.PEAK_PICKER_MIN_INTENSITY_HELP)
            peak_picker_mass_tolerance = c2.number_input(
                label='Peak Picker Mass Tolerance',
                value=qp.peak_picker_mass_tolerance,
                min_value=0.0,
                max_value=1e9,
                help=constants.PEAK_PICKER_MASS_TOLERANCE_HELP)

        with st.expander('Spectra Options'):
            c1, c2 = st.columns(2)
            min_mz = c1.number_input(
                label='Min m/z',
                value=qp.min_mz,
                min_value=0.0,
                max_value=1e9,
                help=constants.MIN_MZ_HELP
            )
            max_mz = c2.number_input(
                label='Max m/z',
                value=qp.max_mz,
                min_value=0.0,
                max_value=1e9,
                help=constants.MAX_MZ_HELP
            )

            c1, c2 = st.columns(2)
            min_intensity_type = c1.selectbox(
                label='Min Intensity Type',
                options=constants.VALID_MIN_INTENSITY_TYPES,
                index=constants.VALID_MIN_INTENSITY_TYPES.index(qp.min_intensity_type),
            )

            min_intensity = c2.number_input(
                label='Min Intensity',
                value=qp.min_intensity,
                min_value=0.0,
                max_value=1e9,
                help=constants.MIN_INTENSITY_HELP
            )

            c1, c2 = st.columns(2)
            top_n = c1.number_input(
                label='Top N',
                value=qp.top_n,
                min_value=0,
                help=constants.TOP_N_HELP,
                disabled=True,
            )

            bottom_n = c2.number_input(
                label='Bottom N',
                value=qp.bottom_n,
                min_value=0,
                help=constants.BOTTOM_N_HELP,
                disabled=True,
            )

            compression_algorithm = st.selectbox(
                label='Compression Algorithm',
                options=constants.VALID_COMPRESSION_ALGORITHMS,
                index=constants.VALID_COMPRESSION_ALGORITHMS.index(qp.compression_algorithm),
            )

        debug = st.checkbox(
            label='Debug',
            value=False
        )
    if spectra:
        mzs, ints = [], []
        for line in spectra.split('\n'):
            mz, intensity = line.split(' ')
            mzs.append(float(mz))
            ints.append(float(intensity))

        int_threshold = min_intensity
        if min_intensity_type == 'relative':
            int_threshold = max(ints) * (min_intensity / 100)

        filtered_mzs, filtered_ints = [], []
        for mz, intensity in zip(mzs, ints):

            if intensity <= int_threshold:
                continue

            if mz < min_mz or mz > max_mz:
                continue

            filtered_mzs.append(mz)
            filtered_ints.append(intensity)

        spectra = list(zip(filtered_mzs, filtered_ints))
    else:
        spectra = []

qp_new = QueryParams(
    sequence=sequence,
    mass_type=mass_type,
    fragment_types=fragment_types,
    mass_tolerance_type=mass_tolerance_type,
    mass_tolerance=mass_tolerance,
    peak_assignment=peak_assignment,
    internal_fragments=internal_fragments,
    min_intensity=min_intensity,
    y_axis_scale=y_axis_scale,
    hide_unassigned_peaks=hide_unassigned_peaks,
    spectra=spectra,
    peak_picker=peak_picker,
    peak_picker_min_intensity=peak_picker_min_intensity,
    peak_picker_mass_tolerance=peak_picker_mass_tolerance,
    neutral_losses=neutral_losses,
    custom_losses=custom_losses,
    isotopes=isotopes,
    filter_missing_mono=filter_missing_mono,
    filter_interrupted_iso=filter_interrupted_iso,
    min_mz=min_mz,
    max_mz=max_mz,
    compression_algorithm=compression_algorithm,
    min_intensity_type=min_intensity_type,
    min_charge=min_charge,
    max_charge=max_charge,
    top_n=top_n,
    bottom_n=bottom_n,
    immonium_ions=immonium_ions,
    internal_fragment_types=internal_fragment_types,
)

if mass_tolerance_type == 'th' and mass_tolerance > 1:
    st.error("Mass tolerance in th must be less than or equal to 1, when type is th.")
    st.stop()

# retrieve ion selections
ion_types = []
charges = []
for i in constants.IONS:
    for c in range(min_charge, max_charge + 1):
        label = get_ion_label(i, c)
        if st.session_state.get(label, False):
            ion_types.append(i)
            charges.append(c)

# retrieve ion selections
internal_ion_types = []
internal_charges = []
for i in constants.INTERNAL_IONS:
    for c in range(min_charge, max_charge + 1):
        label = get_ion_label(i, c)
        if st.session_state.get(label, False):
            internal_ion_types.append(i)
            internal_charges.append(c)

# Show Analysis URL with improved aesthetics
url = generate_app_url(qp_new, qp.compression_algorithm, debug=debug)
url_chars = len(url)
c1, c2 = st.columns([7, 3])
# Use an emoji to make the link more noticeable
c1.write(f'##### :link: [Sharable URL]({url})')

# Add an emoji to signify information about URL length
c2.caption(f'Url Length: {url_chars} characters, 📏 {round(url_chars / 1024, 2)} KB')

st.markdown('---')

# Show Sequence Info
st.header(sequence)
c1, c2, c3, c4 = st.columns(4)

c1.metric('Mass', round(pt.mass(sequence), 4))
c2.metric('Length', len(unmodified_sequence))
c3.metric('Peaks', len(spectra))

cols = st.columns(max_charge - min_charge + 1)
for i in range(min_charge, max_charge + 1):
    col = cols[i - min_charge]
    col.metric(f'M/Z +{i}', round(pt.mz(sequence, charge=i), 4))


@st.cache_data
def build_fragments_cached(*args, **kwargs) -> List[pt.Fragment]:
    return pt.fragment(*args, **kwargs)


annotation = pt.parse(sequence)
fragmenter = pt.Fragmenter(annotation, mass_type == 'monoisotopic')

fragments = build_fragments_cached(sequence=annotation,
                       ion_types=['a', 'b', 'c', 'x', 'y', 'z'],
                       charges=list(range(min_charge, max_charge + 1)),
                       monoisotopic=(mass_type == 'monoisotopic'),
                       isotopes=list(range(isotopes + 1)),
                       losses=losses)

ion_labels = {(ion_type, charge) for ion_type, charge in zip(ion_types, charges)}
fragments = [fragment for fragment in fragments if (fragment.ion_type, fragment.charge) in ion_labels]


if internal_fragments:
    for ion_type, charge in zip(internal_ion_types, internal_charges):
        fragments.extend(fragmenter.fragment(
                                                ion_types=ion_type,
                                                charges=charge,
                                                isotopes=list(range(isotopes + 1)),
                                                losses=losses))

if immonium_ions is True:
    fragments.extend(fragmenter.fragment(
                                                ion_types=['i'],
                                                charges=[1],
                                                isotopes=list(range(isotopes + 1)),
                                                losses=losses))
c4.metric('Fragments', len(fragments))

frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

if spectra:

    mzs, ints = zip(*spectra)

    if peak_picker:
        peaks = [(mz, intensity) for mz, intensity in zip(mzs, ints)]
        deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(peaks, averagine=ms_deisotope.peptide,
                                                               scorer=ms_deisotope.MSDeconVFitter(
                                                                   peak_picker_min_intensity,
                                                                   peak_picker_mass_tolerance))
        mzs = [float(peak.mz) for peak in deconvoluted_peaks]
        ints = [float(peak.intensity) for peak in deconvoluted_peaks]

    # Take top n peaks and bottom n peaks
    #top_n_spectra = sorted(zip(mzs, ints), key=lambda x: x[1], reverse=True)[:top_n]
    #bottom_n_spectra = sorted(zip(mzs, ints), key=lambda x: x[1])[:bottom_n]

    # Combine ensuring no duplicates - since it's a list of tuples, we can use a dictionary to remove duplicates
    # efficiently
    #spectra_dict = {mz: intensity for mz, intensity in top_n_spectra + bottom_n_spectra}
    spectra_dict = {mz: intensity for mz, intensity in zip(mzs, ints)}
    spectra = list(spectra_dict.items())
    mzs, ints = zip(*spectra)

    max_spectra_mz = max(mzs)

    # TODO: Add priority to fragment matches, a random isotope match should not be better than a non-isotope match
    fragment_matches = pt.get_fragment_matches(fragments, mzs, ints, mass_tolerance, mass_tolerance_type,
                                               'largest' if peak_assignment == 'most intense' else 'closest')
    fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)

    if filter_missing_mono:
        fragment_matches = pt.filter_missing_mono_isotope(fragment_matches)

    if filter_interrupted_iso:
        fragment_matches = pt.filter_skipped_isotopes(fragment_matches)

    match_cov = pt.get_match_coverage(fragment_matches)

    if len(fragment_matches) == 0:
        st.warning('No matches found, try increasing the mass tolerance or changing the ion types and charges')
        st.stop()

    fragment_matches = {fm.mz: fm for fm in fragment_matches}  # keep the best fragment match for each mz

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
    spectra_df['matched'] = False
    spectra_df['abs_error'] = spectra_df['error'].abs()
    spectra_df['abs_error_ppm'] = spectra_df['error_ppm'].abs()

    match_df = pd.DataFrame(match_data)
    match_df['matched'] = True
    match_df['abs_error'] = match_df['error'].abs()
    match_df['abs_error_ppm'] = match_df['error_ppm'].abs()

    # spectra_df = spectra_df[~spectra_df['mz'].isin(match_df['mz'])]
    spectra_df = pd.concat([spectra_df, match_df])

    if hide_unassigned_peaks:
        spectra_df = spectra_df[spectra_df['matched']]

    spectra_df['ion_color_type'] = spectra_df['ion_type']
    spectra_df.loc[spectra_df['internal'], 'ion_color_type'] = 'i'

    def create_labels(row):
        if row['ion_type'] != '':
            charge_str = '+' * int(row['charge'])
            charge = int(row['charge'])
            ion_type_str = row['ion_type']

            if row['internal']:
                color_label = get_ion_label(row['ion_type'], int(row['charge']))
                # ion_label = f"{get_ion_label(row['ion_type'], int(row['charge']))}{int(row['parent_number'])}i"
                ion_label = row['label']
            else:
                color_label = get_ion_label(row['ion_type'], int(row['charge']))
                # ion_label = f"{get_ion_label(row['ion_type'], int(row['charge']))}{int(row['parent_number'])}"
                ion_label = row['label']

        else:
            color_label = 'unassigned'
            ion_label = 'unassigned'

        return ion_label, color_label


    # List comprehension to create ion and color labels
    labels = [create_labels(row) for _, row in spectra_df.iterrows()]
    spectra_df['ion_label'], spectra_df['color_label'] = zip(*labels)

    spectra_df.loc[spectra_df['ion_type'] == 'I', 'label'] = spectra_df.loc[
        spectra_df['ion_type'] == 'I', 'sequence'].values

    # Assigning colors based on color labels
    spectra_df['color'] = [color_dict[label] for label in spectra_df['color_label']]

    # change internal ions color to magenta
    # spectra_df.loc[spectra_df['internal'], 'color'] = 'magenta'

    cmap = mpl.colormaps.get_cmap('Blues')

    # Create a color map for the intensities
    st.caption('Sequence coverage')
    for ion, charge in zip(ion_types, charges):
        cov_arr = [0] * len(unmodified_sequence)
        tmp_df = spectra_df[(spectra_df['ion_type'] == ion) & (spectra_df['charge'] == charge)]

        if ion in 'abc':
            for num in tmp_df['end'].unique():
                cov_arr[num - 1] = 1
        elif ion in 'xyz':
            for num in tmp_df['start'].unique():
                cov_arr[num] = 1
        else:
            continue

        if len(cov_arr) > 0:
            c = color_dict[get_ion_label(ion, charge)]
            s = coverage_string(cov_arr, unmodified_sequence, c)

            # center text
            ion_span = f'<span style="color:{c}">{ion}<sup>+{charge}</sup></span>'
            st.markdown(f'{ion_span} {s}', unsafe_allow_html=True)

    fig = generate_annonated_spectra_plotly(spectra_df, scale=y_axis_scale)
    st.plotly_chart(fig, use_container_width=True)

    with tempfile.NamedTemporaryFile(delete=False, suffix='.svg') as tmpfile:
        # Save the figure to the temporary file
        fig.write_image(file=tmpfile.name, format="svg", width=1920, height=1080, scale=3.0)

        # Read the content of the temporary file
        tmpfile.seek(0)  # Go to the start of the file
        data = tmpfile.read()

    # Create a download button in Streamlit
    st.download_button(label="Download chart as SVG",
                       data=data,
                       file_name="spectra.svg",
                       mime="image/svg+xml")

    dfs = []
    # combined_data = {'AA': list(unmodified_sequence)}
    combined_data = {'AA': pt.split(sequence)}
    for ion, charge in zip(ion_types, charges):
        data = {'AA': pt.split(sequence)}
        ion_df = frag_df.copy()
        ion_df = ion_df[
            (ion_df['ion_type'] == ion) & (ion_df['charge'] == charge) &
            (ion_df['internal'] == False) & (ion_df['isotope'] == 0) & (ion_df['loss'] == 0)
            ]
        ion_df.sort_values(by=['start'] if ion in 'xyz' else ['end'], inplace=True,
                           ascending=False if ion in 'xyz' else True)

        # keep only a single number
        ion_df.drop_duplicates(subset=['start'] if ion in 'xyz' else ['end'], inplace=True)

        frags = ion_df['mz'].tolist()

        if ion in 'xyz':
            frags = frags[::-1]

        data[ion] = frags

        combined_data[get_ion_label(ion, charge)] = frags

        # Displaying the table
        df = pd.DataFrame(data)
        df['# (abc)'] = list(range(1, len(df) + 1))
        df['# (xyz)'] = list(range(1, len(df) + 1))[::-1]

        # reorder columns so that # is first # +1 is last and AA is in the middle
        df = df[
            ['AA'] + ['# (abc)'] + [col for col in df.columns if col not in ['AA', '# (abc)', '# (xyz)']] + ['# (xyz)']]
        dfs.append(df)

    combined_df = pd.DataFrame(combined_data)
    # sort columns based on alphabetical order
    combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

    styled_dfs = []


    def color_by_ion_type(col):
        ion_type = col.name[-1]
        color = color_dict.get(ion_type, 'grey')  # get color or default to grey if not found
        return ['color: %s' % color] * len(col)


    for df in dfs:
        styled_df = df.style.apply(color_by_ion_type)

        # Set table styles with increased horizontal padding for more space between columns,
        # centered text, and no borders
        styles = [
            dict(selector="td", props=[("padding", "2px 2px"), ("text-align", "center"), ("border", "none")]),
            dict(selector="th", props=[("padding", "2px 2px"), ("text-align", "center"), ("border", "none")])
        ]
        styled_df = styled_df.set_table_styles(styles)
        styled_dfs.append(styled_df)


    def highlight_cells(data):
        # Initialize empty DataFrame with same index and columns as original
        styled = pd.DataFrame('', index=data.index, columns=data.columns)

        # Iterate over cells and update `styled` based on cell position
        for row in data.index:
            for col in data.columns:
                if col == 'AA' or col == '# (abc)' or col == '# (xyz)':
                    styled.loc[
                        row, col] = f'background-color: gainsboro; color: black; text-align: center; font-weight: bold;'
                    continue

                ion = col[-1]
                if ion in 'abc':
                    ion_number = row + 1
                else:
                    ion_number = len(unmodified_sequence) - row
                label = col + str(ion_number)
                mz = data.loc[row, col]

                if mz <= min_mz or mz >= max_mz:
                    styled.loc[
                        row, col] = f'background-color: #BEBEBE; color: black; text-align: center; font-weight: bold;'
                else:
                    if label in accepted_normal_ions:
                        styled.loc[
                            row, col] = f'background-color: {color_dict[col]}; color: white; text-align: center; font-weight: bold;'
                    elif label in accepted_internal_ions:
                        styled.loc[
                            row, col] = f'background-color: {color_dict[col]}; color: magenta; text-align: center; font-style: italic; font-weight: bold;'
                    else:
                        styled.loc[row, col] = f'background-color: white; color: black; text-align: center;'

        return styled


    matched_ions = spectra_df[spectra_df['ion_type'] != '']
    accepted_normal_ions = matched_ions[matched_ions['internal'] == False]['ion_label'].tolist()
    accepted_internal_ions = matched_ions[matched_ions['internal'] == True]['ion_label'].tolist()
    accepted_internal_ions = [ion[:-1] for ion in accepted_internal_ions]

    combined_df['# (abc)'] = list(range(1, len(unmodified_sequence) + 1))
    combined_df['# (xyz)'] = list(range(1, len(unmodified_sequence) + 1))[::-1]

    # reorder columns so that # is first # +1 is last and AA is in the middle
    combined_cols = combined_df.columns.tolist()
    combined_cols.remove('# (abc)')
    combined_cols.remove('# (xyz)')
    combined_cols.remove('AA')
    forward_cols = [col for col in combined_cols if 'a' in col or 'b' in col or 'c' in col]
    reverse_cols = [col for col in combined_cols if 'x' in col or 'y' in col or 'z' in col]

    # sort
    forward_cols.sort()
    reverse_cols.sort(reverse=True)

    new_cols = ['# (abc)'] + forward_cols + ['AA'] + reverse_cols + ['# (xyz)']
    combined_df = combined_df[new_cols]
    len_combined_df = len(combined_df)

    combined_df = combined_df.style.format(precision=4).apply(highlight_cells, axis=None)

    st.markdown('---')

    st.subheader('Fragment Ions')

    st.dataframe(combined_df, height=int(35.2 * (len_combined_df + 1)), hide_index=True)

    st.plotly_chart(generate_fragment_plot_ion_type(unmodified_sequence, spectra_df), use_container_width=True)

    st.markdown('---')

    st.subheader('Stats')

    # calculate percentage if intensity accounted for
    total_intensity = spectra_df['intensity'].sum()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(label='Total Intensity', value=round(total_intensity, 1))
    c2.metric(label='Matched Intensity', value=round(match_df['intensity'].sum(), 1))
    c3.metric(label='Unmatched Intensity', value=round(spectra_df['intensity'].sum() - match_df['intensity'].sum(), 1))
    c4.metric(label='Matched Intensity %', value=round(match_df['intensity'].sum() / total_intensity * 100, 2))

    if internal_fragments:
        st.plotly_chart(generate_fragment_plot(unmodified_sequence, spectra_df, True), use_container_width=True)

    st.markdown('---')

    # error_fig = generate_error_histogram(spectra_df, mass_tolerance_type)
    # st.plotly_chart(error_fig, use_container_width=True)

    with st.expander('Fragments'):
        st.dataframe(frag_df)

    with st.expander('Peaks'):
        st.dataframe(spectra_df)


else:
    st.warning('No spectra....')
