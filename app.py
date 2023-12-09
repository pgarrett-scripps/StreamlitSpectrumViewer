# TODO: move plotly to separate file
# TODO: push git / streamlit
# TODO: Update and distribute user script
# TODO: Add loss analysis plot

import tempfile

import pandas as pd
import peptacular.sequence
import streamlit as st
from peptacular.fragment import build_fragments
from peptacular.mass import calculate_mass, calculate_mz
from peptacular.score import compute_fragment_matches

import matplotlib as mpl
import ms_deisotope

import constants
from plot_util import generate_annonated_spectra_plotly, coverage_string, generate_fragment_plot, \
    generate_fragment_plot_ion_type
from query_params import parse_query_params, InvalidQueryParam, generate_app_url, QueryParams

st.set_page_config(page_title="Spectra Viewer", page_icon=":glasses:", layout="wide")

params = st.experimental_get_query_params()
try:
    qp = parse_query_params(params)
except InvalidQueryParam as e:
    st.error(str(e))
    st.stop()

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

with st.sidebar:
    st.title('Spectra Viewer')

    sequence = st.text_input(label='Sequence', value=qp.sequence, help=constants.SEQUENCE_HELP)
    unmodified_sequence = peptacular.sequence.strip_modifications(sequence)

    st.write('Ions')
    fragment_types = set()
    for ion in constants.IONS:
        cols = st.columns(constants.MAX_CHARGE)
        for i, col in enumerate(cols, 1):
            label = '+' * i + ion
            if col.checkbox(label,
                            value=label in qp.fragment_types,
                            key=label):
                fragment_types.add(label)

    internal_fragments = st.checkbox(
        label='Internal Fragments',
        value=qp.internal_fragments,
        help=constants.INTERNAL_FRAGMENTS_HELP)

    st.caption('Neutral Losses')
    cols = st.columns(len(constants.NEUTRAL_LOSSES))
    losses = [0.0]
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
        value=';'.join([str(l) for l in qp.custom_losses]),
        help=constants.CUSTOM_LOSSES_HELP
    )

    if custom_losses:
        custom_losses = [float(x) for x in custom_losses.split(';')]
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

    with st.expander('Plot Options'):

        c1, c2 = st.columns(2)
        min_intensity = c1.number_input(
            label='Min Intensity',
            value=qp.min_intensity,
            min_value=0.0,
            max_value=1e9,
            help=constants.MIN_INTENSITY_HELP
        )
        y_axis_scale = c2.radio(
            label='Y Axis Scale',
            options=['linear', 'log'],
            horizontal=True,
            index=0 if qp.y_axis_scale == 'linear' else 1,
            help=constants.Y_AXIS_SCALE_HELP
        )

        hide_unassigned_peaks = st.checkbox(
            label='Hide Unassigned Peaks',
            value=qp.hide_unassigned_peaks,
            help=constants.HIDE_UNASSIGNED_PEAKS_HELP
        )

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

    spectra = st.text_area(
        label='Spectra',
        value='\n'.join(f'{s[0]} {s[1]}' for s in qp.spectra),
        help=constants.SPECTRA_HELP
    )

    if spectra:
        mzs, ints = [], []

        for line in spectra.split('\n'):
            mz, intensity = line.split(' ')

            mz = float(mz)
            intensity = float(intensity)

            if intensity <= min_intensity:
                continue

            if mz < min_mz or mz > max_mz:
                continue

            mzs.append(mz)
            ints.append(intensity)

        spectra = list(zip(mzs, ints))
    else:
        spectra = []

qp = QueryParams(
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
    neutral_losses=set(neutral_losses),
    custom_losses=set(custom_losses),
    isotopes=isotopes,
    filter_missing_mono=filter_missing_mono,
    filter_interrupted_iso=filter_interrupted_iso,
    min_mz=min_mz,
    max_mz=max_mz
)

if mass_tolerance_type == 'th' and mass_tolerance > 1:
    st.error("Mass tolerance in th must be less than or equal to 1, when type is th.")
    st.stop()

# retrieve ion selections
ion_types = []
charges = []
for ion in constants.IONS:
    for i in range(1, constants.MAX_CHARGE + 1):
        label = '+' * i + ion
        if st.session_state.get(label, False):
            ion_types.append(ion)
            charges.append(i)

# Show Analysis URL
url = generate_app_url(qp)
st.write(f'##### [Analysis URL]({url}) (copy me and send to your friends!)')

# Show Sequence Info
sequence_charge = 2
st.header(sequence)
c1, c2, c3, c4 = st.columns(4)
c1.metric('Mass', round(calculate_mass(sequence), 4))
c2.metric('m/z', round(calculate_mz(sequence, sequence_charge), 4))
c3.metric('Charge', sequence_charge)
c4.metric('Length', len(unmodified_sequence))

fragments = []
for ion, charge in zip(ion_types, charges):
    fragments.extend(build_fragments(sequence=sequence,
                                     ion_types=ion,
                                     charges=charge,
                                     monoisotopic=(mass_type == 'monoisotopic'),
                                     internal=internal_fragments,
                                     isotopes=list(range(isotopes + 1)),
                                     losses=losses))

frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

if spectra:

    mzs, ints = zip(*spectra)

    if peak_picker:
        peaks = [(mz, intensity) for mz, intensity in zip(mzs, ints)]
        deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(peaks, averagine=ms_deisotope.peptide,
                                                               scorer=ms_deisotope.MSDeconVFitter(
                                                                   qp.peak_picker_min_intensity,
                                                                   qp.peak_picker_mass_tolerance))
        mzs = [float(peak.mz) for peak in deconvoluted_peaks]
        ints = [float(peak.intensity) for peak in deconvoluted_peaks]

    max_spectra_mz = max(mzs)

    # TODO: Add priority to fragment matches, a random isotope match should not be better than a non-isotope match
    fragment_matches = compute_fragment_matches(fragments, mzs, ints, mass_tolerance, mass_tolerance_type)
    fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)
    fragment_matches = {fm.mz: fm for fm in fragment_matches}  # keep the best error for each fragment

    match_data = {'sequence': [], 'charge': [], 'ion_type': [], 'number': [], 'internal': [], 'parent_number': [],
                  'monoisotopic': [], 'mz': [], 'intensity': [], 'error': [], 'abs_error': [], 'theo_mz': [],
                  'label': [], 'isotope': [], 'loss': [], 'error_ppm': [], 'abs_error_ppm': [], 'fragment_mz': [],
                  'start': [], 'end': []}
    data = {'sequence': [], 'charge': [], 'ion_type': [], 'number': [], 'internal': [], 'parent_number': [],
            'monoisotopic': [], 'mz': [], 'intensity': [], 'error': [], 'abs_error': [], 'theo_mz': [],
            'label': [], 'isotope': [], 'loss': [], 'error_ppm': [], 'abs_error_ppm': [], 'fragment_mz': [],
            'start': [], 'end': []}

    for mz, i in zip(mzs, ints):
        fm = fragment_matches.get(mz, None)

        if fm:
            match_data['sequence'].append(fm.fragment.sequence)
            match_data['charge'].append(fm.fragment.charge)
            match_data['ion_type'].append(fm.fragment.ion_type)
            match_data['number'].append(fm.fragment.number)
            match_data['internal'].append(fm.fragment.internal)
            match_data['parent_number'].append(fm.fragment.parent_number)
            match_data['monoisotopic'].append(fm.fragment.monoisotopic)
            match_data['mz'].append(mz)
            match_data['intensity'].append(i)
            match_data['error'].append(fm.error)
            match_data['abs_error'].append(abs(fm.error))
            match_data['error_ppm'].append(fm.error_ppm)
            match_data['abs_error_ppm'].append(abs(fm.error_ppm))
            match_data['theo_mz'].append(fm.fragment.mz)
            match_data['label'].append(fm.fragment.label)
            match_data['isotope'].append(fm.fragment.isotope)
            match_data['loss'].append(fm.fragment.loss)
            match_data['fragment_mz'].append(fm.fragment.mz)
            match_data['start'].append(fm.fragment.start)
            match_data['end'].append(fm.fragment.end)

        data['sequence'].append('')
        data['charge'].append(0)
        data['ion_type'].append('')
        data['number'].append(0)
        data['internal'].append(False)
        data['parent_number'].append(0)
        data['monoisotopic'].append(True)
        data['mz'].append(mz)
        data['intensity'].append(i)
        data['error'].append(0)
        data['abs_error'].append(0)
        data['error_ppm'].append(0)
        data['abs_error_ppm'].append(0)
        data['theo_mz'].append(0)
        data['label'].append('')
        data['isotope'].append(0)
        data['loss'].append(0.0)
        data['fragment_mz'].append(0.0)
        data['start'].append(0)
        data['end'].append(0)

    spectra_df = pd.DataFrame(data)
    spectra_df['matched'] = False
    match_df = pd.DataFrame(match_data)
    match_df['matched'] = True

    # for keep only the lowest abs_error for ion_type, charge, num
    if peak_assignment == 'most intense':
        match_df.sort_values(by='intensity', inplace=True, ascending=False)
        match_df.drop_duplicates(subset=['theo_mz'], inplace=True)

    else:
        match_df.sort_values(by='abs_error', inplace=True)
        match_df.drop_duplicates(subset=['theo_mz'], inplace=True)

    if filter_missing_mono:
        # remove peaks that skip isotopes
        mono_labels = set(match_df[match_df['isotope'] == 0]['label'].unique())
        labels_to_remove = set()
        for label in match_df[match_df['isotope'] != 0]['label'].unique():
            mono_label = label.replace('*', '')
            if mono_label not in mono_labels:
                labels_to_remove.add(label)

        match_df = match_df[~match_df['label'].isin(labels_to_remove)]

    if filter_interrupted_iso:
        # remove peaks any peaks after a missing isotope
        labels_to_remove = set()

        # get min isotope label for each mono label
        labels = set(match_df['label'].unique())
        mono_iso_labels = set([label.replace('*', '') for label in labels])
        min_iso_labels = set()
        for mono_label in mono_iso_labels:
            for i in range(0, isotopes + 1):
                label = mono_label + '*' * i
                if label in labels:
                    min_iso_labels.add(label)
                    break

        for min_label in min_iso_labels:
            break_flag = False
            iso_cnt = min_label.count('*')
            for i in range(iso_cnt, isotopes + 1):
                label = min_label + '*' * i
                if label not in labels:
                    break_flag = True
                if break_flag:
                    labels_to_remove.add(label)
        match_df = match_df[~match_df['label'].isin(labels_to_remove)]

    spectra_df = spectra_df[~spectra_df['mz'].isin(match_df['mz'])]
    spectra_df = pd.concat([spectra_df, match_df])

    if hide_unassigned_peaks:
        spectra_df = spectra_df[spectra_df['matched']]

    spectra_df['ion_color_type'] = spectra_df['ion_type']
    spectra_df.loc[spectra_df['internal'], 'ion_color_type'] = 'i'


    def create_labels(row):
        if row['ion_type'] != '':
            charge_str = '+' * int(row['charge'])
            ion_type_str = row['ion_type']
            parent_number_str = str(int(row['parent_number']))
            internal_str = 'i' if row['internal'] else ''

            if row['internal']:
                color_label = f"{charge_str}{internal_str}"
            else:
                color_label = f"{charge_str}{ion_type_str}"

            ion_label = f"{charge_str}{ion_type_str}{parent_number_str}{internal_str}"
        else:
            color_label = 'unassigned'
            ion_label = 'unassigned'

        return ion_label, color_label


    # List comprehension to create ion and color labels
    labels = [create_labels(row) for _, row in spectra_df.iterrows()]
    spectra_df['ion_label'], spectra_df['color_label'] = zip(*labels)

    # Assigning colors based on color labels
    spectra_df['color'] = [constants.COLOR_DICT[label] for label in spectra_df['color_label']]

    # change internal ions color to magenta
    #spectra_df.loc[spectra_df['internal'], 'color'] = 'magenta'

    cmap = mpl.colormaps.get_cmap('Blues')

    st.caption('Sequence coverage')
    for ion, charge in zip(ion_types, charges):
        cov_arr = [0] * len(unmodified_sequence)
        tmp_df = spectra_df[(spectra_df['ion_type'] == ion) & (spectra_df['charge'] == charge)]
        nums = tmp_df['parent_number'].unique()

        if ion in 'abc':
            for num in nums:
                cov_arr[num - 1] = 1
        else:
            for num in nums:
                cov_arr[len(unmodified_sequence) - (num - 1) - 1] = 1

        if len(cov_arr) > 0:
            c = constants.COLOR_DICT['+' * charge + ion]
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
    combined_data = {'AA': list(unmodified_sequence)}
    for ion, charge in zip(ion_types, charges):
        data = {'AA': list(unmodified_sequence)}
        ion_df = frag_df[
            (frag_df['ion_type'] == ion) & (frag_df['charge'] == charge) &
            (frag_df['internal'] == False)]
        ion_df.sort_values(by=['number'], inplace=True)

        # keep only a single number
        ion_df.drop_duplicates(subset=['number'], inplace=True)

        frags = ion_df['mz'].tolist()

        if ion in 'xyz':
            frags = frags[::-1]

        data[ion] = frags

        combined_data['+' * charge + ion] = frags

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
        color = constants.COLOR_DICT.get(ion_type, 'grey')  # get color or default to grey if not found
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
                if label in accepted_normal_ions:
                    styled.loc[
                        row, col] = f'background-color: {constants.COLOR_DICT[col]}; color: white; text-align: center; font-weight: bold;'
                elif label in accepted_internal_ions:
                    styled.loc[
                        row, col] = f'background-color: {constants.COLOR_DICT[col]}; color: magenta; text-align: center; font-style: italic; font-weight: bold;'
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

    # calculate percerntage if intensity accounted for
    total_intensity = spectra_df['intensity'].sum()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(label='Total Intensity', value=round(total_intensity, 1))
    c2.metric(label='Matched Intensity', value=round(match_df['intensity'].sum(), 1))
    c3.metric(label='Unmatched Intensity', value=round(spectra_df['intensity'].sum() - match_df['intensity'].sum(), 1))
    c4.metric(label='Matched Intensity %', value=round(match_df['intensity'].sum() / total_intensity * 100, 2))

    if internal_fragments:
        st.plotly_chart(generate_fragment_plot(unmodified_sequence, spectra_df, True), use_container_width=True)

    st.markdown('---')

    st.subheader('Isotopes')

    iso_data = []
    for i, row in spectra_df.iterrows():
        mz = row['mz']
        intensity = row['intensity']
        charge = row['charge']

        if charge == 0:
            continue

        isotopic_pattern = ms_deisotope.peptide.isotopic_cluster(mz, charge)
        for peak in isotopic_pattern:
            iso_data.append([i, peak.mz, peak.intensity])

with st.expander('Fragments'):
    st.dataframe(frag_df)

with st.expander('Peaks'):
    st.dataframe(spectra_df)
