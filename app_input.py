from functools import cached_property

import streamlit as st
import streamlit_permalink as stp
import peptacular as pt
from dataclasses import dataclass
from typing import List, Tuple

import constants
from color_util import get_color_dict
from msms_compression import SpectrumCompressorUrl
from msdecon.deconvolution import deconvolute

@dataclass
class SpectraInputs:
    """Dataclass to store all spectra viewer inputs."""
    # Sequence parameters
    sequence: str
    
    # Mass tolerance parameters
    mass_tolerance_type: str
    mass_tolerance: float
    
    # Charge parameters
    min_charge: int
    max_charge: int
    
    # Spectra data
    spectra_text: str
    
    # Fragment parameters
    fragment_types: List[str]
    immonium_ions: bool
    mass_type: str
    peak_assignment: str
    
    # Isotope parameters
    num_isotopes: int
    filter_missing_mono: bool
    filter_interrupted_iso: bool
    
    # Plot parameters
    y_axis_scale: str
    hide_unassigned_peaks: bool
    min_mz: float
    max_mz: float
    min_intensity_type: str
    min_intensity: float
    max_intensity_type: str
    max_intensity: float
    line_width: float
    text_size: float
    marker_size: float
    axis_text_size: float
    title_text_size: float
    tick_text_size: float
    fig_width: float
    fig_height: float
    hide_error_percentile_labels: bool
    bold_labels: bool
    color_dict: dict[str, str]
    
    # Neutral loss parameters
    h2o_loss: bool
    nh3_loss: bool
    h3po4_loss: bool
    custom_loss_str: str

    #deconvolution parameters
    deconvolute: bool
    deconvolute_error_type: str
    deconvolute_error: float

    stateful: bool = True
    @property
    def custom_losses(self) -> dict[str, float]:

        if self.custom_loss_str:
            return {
                nl.split(":")[0]: float(nl.split(":")[1])
                for nl in self.custom_loss_str.split(";")
            }
        else:
            return {}
    
    @property
    def neutral_losses(self) -> dict[str, float]:
        losses = {}
        if self.h2o_loss:
            losses['[STED]'] = -18.01056
        if self.nh3_loss:
            losses['[RKNQ]'] = -17.02655
        if self.h3po4_loss:
            losses['[ST]'] = -97.9769

        losses.update(self.custom_losses)

        return losses
    
    @property
    def losses(self) -> list[tuple[str, float]]:
        return list(self.neutral_losses.items()) + list(self.custom_losses.items())
    
    @property
    def ion_types(self) -> list[str]:
        return [
            f"{f}{c}"
            for f in self.fragment_types
            for c in range(self.min_charge, self.max_charge + 1)
        ]

    @cached_property
    def spectra(self) -> list[tuple[float, float]]:

        spectra = parse_sequence(self.spectra_text)

        # filter
        if self.min_intensity_type == "relative":
            max_intensity = max([intensity for _, intensity in spectra])
            spectra = [
                (mz, intensity)
                for mz, intensity in spectra
                if intensity >= self.min_intensity / 100 * max_intensity
            ]

        if self.min_intensity_type == "absolute":
            spectra = [
                (mz, intensity)
                for mz, intensity in spectra
                if intensity >= self.min_intensity
            ]

        if self.max_intensity_type == "relative":
            max_intensity = max([intensity for _, intensity in spectra])
            spectra = [
                (mz, intensity)
                for mz, intensity in spectra
                if intensity <= self.max_intensity / 100 * max_intensity
            ]

        if self.max_intensity_type == "absolute":
            spectra = [
                (mz, intensity)
                for mz, intensity in spectra
                if intensity <= self.max_intensity
            ]

        if self.min_mz:
            spectra = [
                (mz, intensity)
                for mz, intensity in spectra
                if mz >= self.min_mz
            ]

        if self.max_mz:
            spectra = [
                (mz, intensity)
                for mz, intensity in spectra
                if mz <= self.max_mz
            ]

        if self.deconvolute:
            peaks = deconvolute(spectra,
                                  tolerance=self.deconvolute_error,
                                  tolerance_type=self.deconvolute_error_type,
                                  charge_range=(self.min_charge, self.max_charge))

            spectra = [(p.base_peak.mz, p.total_intensity) for p in peaks]

        return spectra

    @cached_property
    def min_spectra_mz(self):
        return min([mz for mz, _ in self.spectra])

    @cached_property
    def max_spectra_mz(self):
        return max([mz for mz, _ in self.spectra])
    
    @cached_property
    def min_spectra_intensity(self):
        return min([intensity for _, intensity in self.spectra])
    
    @cached_property
    def max_spectra_intensity(self):
        return max([intensity for _, intensity in self.spectra])
    
    @cached_property
    def mz_values(self) -> list[float]:
        return [mz for mz, _ in self.spectra]
    
    @cached_property
    def intensity_values(self) -> list[float]:
        return [intensity for _, intensity in self.spectra]

    @cached_property
    def mz_int_values(self) -> (list[float], list[float]):
        return self.mz_values, self.intensity_values
    
    
    def get_color(self, ion: str, charge: int) -> str:
        """Get color for a specific ion and charge."""

        if ion is None or charge is None:
            return self.color_dict['unassigned']

        color_key = f'{"+" * charge}{ion}'
        if color_key in self.color_dict:
            return self.color_dict[color_key]
        else:
            return self.color_dict['unassigned']

    
    @property
    def charges(self) -> list[int]:
        return [
            c
            for c in range(self.min_charge, self.max_charge + 1)
        ]
    

    @property
    def unmodified_sequence(self) -> str:
        return pt.strip_mods(self.sequence)

    @property
    def is_monoisotopic(self) -> bool:
        return self.mass_type == "monoisotopic"

    @property
    def isotopes(self) -> list[int]:
        return list(range(self.num_isotopes + 1))

    @property
    def peak_assignment_type(self) -> str:
        return "largest" if self.peak_assignment == "most intense" else "closest"


def get_ion_label(i: str, c: int) -> str:
    """Get ion label with charge."""
    return "+" * c + i


def parse_sequence(input_str: str) -> List[Tuple[float, float]]:
    """Parse sequence string into list of mz and intensity tuples."""

    if not input_str:
        return []

    mz_values, intensity_values = [], []
    lines = input_str.split("\n")
    for line in lines:
        parts = line.split(" ")
        if len(parts) == 2:
            try:
                mz_values.append(float(parts[0]))
                intensity_values.append(float(parts[1]))
            except ValueError:
                raise ValueError(f"Invalid input line: {parts}")
        else:
            raise ValueError(f"Invalid input format: {input_str}")
    return list(zip(mz_values, intensity_values))


def serialize_sequence(sequence: List[Tuple[float, float]]) -> str:
    """Convert list of mz and intensity tuples to string format."""
    return "\n".join(f"{round(s[0], 5)} {round(s[1], 2)}" for s in sequence)

@st.cache_data
def compress_spectra(input_str: str) -> str:
    """Compress spectra string for URL encoding."""
    try:
        spectra = parse_sequence(input_str)
        mzs, ints = [], []
        if spectra:
            mzs, ints = zip(*spectra)
        value = SpectrumCompressorUrl.compress(mzs, ints)
        return value
    except ValueError as e:
        st.error(f"Error compressing spectra: {e}")
        return SpectrumCompressorUrl.compress([], [])


@st.cache_data
def decompress_spectra(input_str: str) -> str:
    """Decompress spectra string from URL encoding."""
    if input_str.startswith("RAW:"):
        mzs, ints = zip(
            *[
                (float(elem.split(":")[0]), float(elem.split(":")[1]))
                for elem in input_str[4:].split(";")
            ]
        )
        return serialize_sequence(list(zip(mzs, ints)))

    try:
        mzs, ints = SpectrumCompressorUrl.decompress(input_str)
        spectra = list(zip(mzs, ints))
        return serialize_sequence(spectra)
    except ValueError as e:
        st.error(f"Error decompressing spectra: {e}")
        raise ValueError(f"Error decompressing spectra: {e}") from e


def get_all_inputs(stateful: bool) -> SpectraInputs:
    """Get all inputs from the Streamlit UI and return as a SpectraInputs dataclass."""

    # Inject custom CSS to set the width of the sidebar
    st.markdown(
        """
        <style>
            section[data-testid="stSidebar"] {
                width: 600px !important; # Set the width to your desired value
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    input_tab, frag_tab, iso_tab, loss_tab, deconv_tab, plot_tab, vis_tab = st.tabs(
        ["Input", "Fragment", "Isotope", "Loss", "Deconv", "Spectra", "Plot"])

    with input_tab:
        # Sequence input
        sequence = stp.text_input(
            label="Sequence (Proforma2.0 Notation)",
            value=constants.DEFAULT_SEQUENCE,
            help=constants.SEQUENCE_HELP,
            key="peptide_sequence",
            stateful=stateful,
        )

        spectra_text = stp.text_area(
            label="Spectra (mz Intensity)",
            value=serialize_sequence(constants.DEFAULT_SPECTRA),
            help=constants.SPECTRA_HELP,
            height=150,
            key="spectra",
            compress=True,
            compressor=compress_spectra,
            decompressor=decompress_spectra,
            stateful=stateful,
        )

    with frag_tab:

        # Fragment ion selection
        fragment_types = stp.pills(
            "Fragment Ions",
            selection_mode="multi",
            options=constants.FRAGMENT_TYPES + ['immonium'],
            default=constants.DEFAULT_FRAGMENT_TYPES,
            help="Select fragment ion types to display",
            key="fragment_types",
            stateful=stateful,
        )

        immonium_ions = 'immonium' in fragment_types

        # drop immonium
        fragment_types = [f for f in fragment_types if f != 'immonium']

        # Mass tolerance settings
        c1, c2 = st.columns(2)
        with c1:
            mass_tolerance_type = stp.selectbox(
                label="Mass Tolerance Type",
                options=constants.MASS_TOLERANCE_TYPES,
                index=constants.MASS_TOLERANCE_TYPES.index(
                    constants.DEFAULT_MASS_TOLERANCE_TYPE
                ),
                help=constants.MASS_TOLERANCE_TYPE_HELP,
                key="mass_tolerance_type",
                stateful=stateful,
            )

        with c2:
            mass_tolerance = stp.number_input(
                label="Mass Tolerance",
                value=constants.DEFAULT_PPM_MASS_TOLERANCE,
                help=constants.MASS_TOLERANCE_HELP,
                key="mass_error",
                stateful=stateful,
            )

        # Charge range
        c1, c2 = st.columns(2)

        with c1:
            min_charge = stp.number_input(
                label="Min Charge",
                min_value=0,
                max_value=100,
                value=1,
                key="min_charge",
                #on_change=on_min_change,
                stateful=stateful,
            )

        with c2:
            max_charge = stp.number_input(
                label="Max Charge",
                min_value=0,
                max_value=100,
                value=2,
                key="max_charge",
                #on_change=on_max_change,
                stateful=stateful,
            )

        # Mass type and peak assignment
        c1, c2 = st.columns(2)
        with c1:
            mass_type = stp.selectbox(
                label="Mass Type",
                options=constants.VALID_MASS_TYPES,
                index=constants.VALID_MASS_TYPES.index(constants.DEFAULT_MASS_TYPE),
                help=constants.MASS_TYPE_HELP,
                key="mass_type",
                stateful=stateful,
            )

        with c2:
            peak_assignment = stp.selectbox(
                label="Peak Assignment",
                options=constants.PEAK_ASSIGNMENTS,
                index=constants.PEAK_ASSIGNMENTS.index(constants.DEFAULT_PEAK_ASSIGNMENT),
                help=constants.PEAK_ASSIGNMENT_HELP,
                key="peak_assignment",
                stateful=stateful,
            )

    # Isotope settings
    with iso_tab:
        num_isotopes = stp.number_input(
            label="Isotopes",
            value=0,
            min_value=0,
            max_value=5,
            help=constants.ISOTOPES_HELP,
            stateful=stateful,
        )

        c1, c2 = st.columns(2)
        with c1:
            filter_missing_mono = stp.checkbox(
                label="Filter missing mono peaks",
                value=False,
                help=constants.FILTER_MISSING_MONO_HELP,
                key="filter_missing_mono",
                stateful=stateful,
            )

        with c2:
            filter_interrupted_iso = stp.checkbox(
                label="Filter interrupted isotopes",
                value=False,
                help=constants.FILTER_INTERRUPTED_ISO_HELP,
                key="filter_interrupted_iso",
                stateful=stateful,
            )

    # Plot options
    with plot_tab:
        c1, c2 = st.columns(2)
        with c1:
            y_axis_scale = stp.radio(
                label="Y Axis Scale",
                options=["linear", "log"],
                horizontal=True,
                index=0,
                help=constants.Y_AXIS_SCALE_HELP,
                stateful=stateful,
            )

        with c2:
            hide_unassigned_peaks = stp.checkbox(
                label="Hide Unassigned Peaks",
                value=False,
                help=constants.HIDE_UNASSIGNED_PEAKS_HELP,
                key="hide_unassigned_peaks",
                stateful=stateful,
            )

        c1, c2 = st.columns(2)
        with c1:
            min_mz = stp.number_input(
                label="Min m/z",
                value=0.0,
                min_value=0.0,
                max_value=1_000_000.0,
                help=constants.MIN_MZ_HELP,
                key="min_mz",
                stateful=stateful,
            )

            min_intensity_type = stp.selectbox(
                label="Min Intensity Type",
                options=constants.VALID_MIN_INTENSITY_TYPES,
                index=constants.VALID_MIN_INTENSITY_TYPES.index(
                    constants.DEFAULT_MIN_INTENSITY_TYPE
                ),
                key="min_intensity_type",
                stateful=stateful,
            )

            max_intensity_type = stp.selectbox(
                label="Max Intensity Type",
                options=constants.VALID_MIN_INTENSITY_TYPES,
                index=constants.VALID_MIN_INTENSITY_TYPES.index(
                    constants.DEFAULT_MIN_INTENSITY_TYPE
                ),
                key="max_intensity_type",
                stateful=stateful,
            )

        with c2:
            max_mz = stp.number_input(
                label="Max m/z",
                value=1_000_000.0,
                min_value=0.0,
                max_value=1e9,
                help=constants.MAX_MZ_HELP,
                key="max_mz",
                stateful=stateful,
            )

            min_intensity = stp.number_input(
                label="Min Intensity",
                value=0.0,
                min_value=0.0,
                max_value=1e9,
                help=constants.MIN_INTENSITY_HELP,
                key="min_intensity",
                stateful=stateful,
            )
            max_intensity = stp.number_input(
                label="Max Intensity",
                value=1e9,
                min_value=0.0,
                max_value=1e9,
                help=constants.MIN_INTENSITY_HELP,
                key="max_intensity",
                stateful=stateful,
            )

    with vis_tab:

        c1, c2 = st.columns(2)
        with c1:
            line_width = stp.number_input(
                label="Line Width",
                value=2.0,
                min_value=0.0,
                max_value=5.0,
                step=0.25,
                help=constants.LINE_WIDTH_HELP,
                key="line_width",
                stateful=stateful,
            )


            text_size = stp.number_input(
                label="Text Size",
                value=20.0,
                min_value=8.0,
                max_value=30.0,
                step=1.0,
                help=constants.TEXT_SIZE_HELP,
                key="text_size",
                stateful=stateful,
            )

            tick_text_size = stp.number_input(
                label="Tick Text Size",
                value=15.0,
                min_value=8.0,
                max_value=30.0,
                step=1.0,
                key="tick_text_size",
                stateful=stateful,
            )
            
            fig_height = stp.number_input(
                label="Figure Height",
                value=600,
                min_value=100,
                max_value=5000,
                step=100,
                key="fig_height",
                stateful=stateful,
            )

            
            bold_labels = stp.checkbox(
                label="Bold Labels",
                value=False,
                key="bold_labels",
                stateful=stateful,
            )
        with c2:
            marker_size = stp.number_input(
                label="Marker Size",
                value=6.0,
                min_value=1.0,
                max_value=30.0,
                step=1.0,
                help=constants.MARKER_SIZE_HELP,
                key="marker_size",
                stateful=stateful,
            )

            axis_text_size = stp.number_input(
                label="Axis Text Size",
                value=15.0,
                min_value=8.0,
                max_value=30.0,
                step=1.0,
                key="axis_text_size",
                stateful=stateful,
            )
            title_text_size = stp.number_input(
                label="Title Text Size",
                value=20.0,
                min_value=8.0,
                max_value=30.0,
                step=1.0,
                key="title_text_size",
                stateful=stateful,
            )


            fig_width = stp.number_input(
                label="Figure Width",
                value=800,
                min_value=100,
                max_value=5000,
                step=100,
                key="fig_width",
                stateful=stateful,
            )


            hide_error_percentile_labels = stp.checkbox(
                label="Hide Error Percentile Labels",
                value=False,
                key="hide_error_percentile_labels",
                stateful=stateful,
            )

        _default_color_dict = get_color_dict(min_charge, max_charge)

        color_dict = {'unassigned': _default_color_dict['unassigned']}
        for ion_type in fragment_types:
            for charge in range(min_charge, max_charge + 1):
                key = f"{'+'*charge}{ion_type}"
                color_dict[key] = _default_color_dict[key]

        with st.expander('Colors'):
            cols = st.columns([1, 1, 1])
            for i, key in enumerate(color_dict.keys()):
                with cols[i % len(cols)]:
                    color_dict[key] = stp.color_picker(f"{key}", value=color_dict[key], stateful=stateful)


    # Neutral losses
    with loss_tab:

        loss_pills = stp.pills('Neutral Losses',
                            selection_mode='multi',
                            options=['H2O', 'NH3', 'H3PO4'],
                            default=[],
                            help='Select neutral losses to display',
                            key='losses',
                            stateful=stateful)

        h2o_loss = 'H2O' in loss_pills
        nh3_loss = 'NH3' in loss_pills
        h3po4_loss = 'H3PO4' in loss_pills
        
        # Handle custom losses
        custom_loss_str = stp.text_input(
            label="Custom Losses",
            value="",
            help=constants.CUSTOM_LOSSES_HELP,
            placeholder="[STED]:-18.01056;[RKNQ]:-17.02655;[ST]:-97.9769",
            key="custom_losses",
            stateful=stateful,
        )

    # Deconvolution options
    with deconv_tab:
        deconvolute = stp.checkbox(label='Deconvolute',
                        value=False,
                        key='deconvolute',
                        stateful=stateful)

        c1, c2 = st.columns(2)
        with c1:
            deconvolute_error_type = stp.selectbox(
                label="Error Type",
                options=['da', 'ppm'],
                index=1,
                key="deconvolute_error_type",
                stateful=stateful,
            )

        with c2:
            deconvolute_error = stp.number_input(
                label="Error",
                value=10.0,
                key="deconvolute_error",
                stateful=stateful,
            )


    # Create and return the dataclass instance
    return SpectraInputs(
        sequence=sequence,
        mass_tolerance_type=mass_tolerance_type,
        mass_tolerance=mass_tolerance,
        spectra_text=spectra_text,
        min_charge=min_charge,
        max_charge=max_charge,
        fragment_types=fragment_types,
        immonium_ions=immonium_ions,
        mass_type=mass_type,
        peak_assignment=peak_assignment,
        num_isotopes=num_isotopes,
        filter_missing_mono=filter_missing_mono,
        filter_interrupted_iso=filter_interrupted_iso,
        y_axis_scale=y_axis_scale,
        hide_unassigned_peaks=hide_unassigned_peaks,
        min_mz=min_mz,
        max_mz=max_mz,
        min_intensity_type=min_intensity_type,
        min_intensity=min_intensity,
        max_intensity_type=max_intensity_type,
        max_intensity=max_intensity,
        h2o_loss=h2o_loss,
        nh3_loss=nh3_loss,
        h3po4_loss=h3po4_loss,
        custom_loss_str=custom_loss_str,
        deconvolute=deconvolute,
        deconvolute_error_type=deconvolute_error_type,
        deconvolute_error=deconvolute_error,
        stateful=stateful,
        text_size=text_size,
        line_width=line_width,
        marker_size=marker_size,
        axis_text_size=axis_text_size,
        title_text_size=title_text_size,
        tick_text_size=tick_text_size,
        fig_width=fig_width,
        fig_height=fig_height,
        hide_error_percentile_labels=hide_error_percentile_labels,
        bold_labels=bold_labels,
        color_dict=color_dict
    )

