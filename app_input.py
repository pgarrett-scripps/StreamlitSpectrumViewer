from functools import cache, lru_cache

import streamlit as st
import streamlit_permalink as stp
import peptacular as pt
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Dict, Any, Optional

import constants
from color_util import get_color_dict
from msms_compression import SpectrumCompressorUrl


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
    
    # Neutral loss parameters
    h2o_loss: bool
    nh3_loss: bool
    h3po4_loss: bool
    custom_loss_str: str

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

    @property
    def spectra(self) -> list[tuple[float, float]]:
        return parse_sequence(self.spectra_text)
    
    @property
    def mz_values(self) -> list[float]:
        return [mz for mz, _ in self.spectra]
    
    @property
    def intensity_values(self) -> list[float]:
        return [intensity for _, intensity in self.spectra]
    
    @property
    def filtered_spectra(self) -> list[tuple[float, float]]:
        return filter_spectra(
            self.spectra_text, 
            self.min_intensity, 
            self.min_intensity_type, 
            self.min_mz, 
            self.max_mz)
    
    @property
    def filtered_mz_values(self) -> list[float]:
        return [mz for mz, _ in self.filtered_spectra]
    
    @property
    def filtered_intensity_values(self) -> list[float]:
        return [intensity for _, intensity in self.filtered_spectra]
    
    @property
    def color_dict(self) -> dict[str, str]:
        return get_color_dict(self.min_charge, self.max_charge)
    
    @property
    def ion_types(self) -> list[str]:
        return [
            f"{f}{c}"
            for f in self.fragment_types
            for c in range(self.min_charge, self.max_charge + 1)
        ]
    
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


def get_ion_label(i: str, c: int) -> str:
    """Get ion label with charge."""
    return "+" * c + i


def parse_sequence(input_str: str) -> List[Tuple[float, float]]:
    """Parse sequence string into list of mz and intensity tuples."""
    mz_values, intensity_values = [], []
    lines = input_str.split("\n")
    for line in lines:
        parts = line.split(" ")
        if len(parts) == 2:
            try:
                mz_values.append(float(parts[0]))
                intensity_values.append(float(parts[1]))
            except ValueError:
                raise ValueError("Invalid input format")
        else:
            raise ValueError("Invalid input format")
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
        raise e


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
        st.error(input_str)
        st.error(f"Error decompressing spectra: {e}")
        raise e


def filter_spectra(spectra, min_intensity, min_intensity_type, min_mz, max_mz):
    """Filter spectra based on intensity and m/z thresholds."""
    if not spectra:
        return []
    
    mzs, ints = [], []
    for line in spectra.split("\n"):
        mz, intensity = line.split(" ")
        mzs.append(float(mz))
        ints.append(float(intensity))

    int_threshold = min_intensity
    if min_intensity_type == "relative":
        int_threshold = max(ints) * (min_intensity / 100)

    filtered_mzs, filtered_ints = [], []
    for mz, intensity in zip(mzs, ints):
        if intensity <= int_threshold:
            continue
        if mz < min_mz or mz > max_mz:
            continue
        filtered_mzs.append(mz)
        filtered_ints.append(intensity)

    return list(zip(filtered_mzs, filtered_ints))


def get_all_inputs(stateful: bool) -> SpectraInputs:
    """Get all inputs from the Streamlit UI and return as a SpectraInputs dataclass."""

    # Sequence input
    sequence = stp.text_input(
        label="Sequence (Proforma2.0 Notation)",
        value=constants.DEFAULT_SEQUENCE,
        help=constants.SEQUENCE_HELP,
        key="peptide_sequence",
        stateful=stateful,
    )

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

    # Fragment ion selection
    fragment_types = stp.pills(
        "Fragment Ions",
        selection_mode="multi",
        options=constants.FRAGMENT_TYPES,
        default=constants.DEFAULT_FRAGMENT_TYPES,
        help="Select fragment ion types to display",
        key="fragment_types",
        stateful=stateful,
    )

    # Immonium ions
    immonium_ions = stp.checkbox(
        label="Immonium Ions",
        value=constants.DEFAULT_IMMONIUM_IONS,
        help=constants.IMMONIUM_IONS_HELP,
        key="immonium_ions",
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
    with st.expander("Isotopes"):
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
    with st.expander("Plot Options"):
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

    # Neutral losses
    with st.expander("Neutral Losses"):

        # Water loss
        h2o_loss = stp.checkbox(label='H2O', 
                        value=False, 
                        key='H2O', 
                        stateful=stateful)
        

        nh3_loss = stp.checkbox(label='NH3', 
                        value=False, 
                        key='NH3', 
                        stateful=stateful)
        
        h3po4_loss = stp.checkbox(label='H3PO4', 
                        value=False, 
                        key='H3PO4', 
                        stateful=stateful)         
        
        # Handle custom losses
        custom_loss_str = stp.text_input(
            label="Custom Losses",
            value="",
            help=constants.CUSTOM_LOSSES_HELP,
            placeholder="[STED]:-18.01056;[RKNQ]:-17.02655;[ST]:-97.9769",
            key="custom_losses",
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
        h2o_loss=h2o_loss,
        nh3_loss=nh3_loss,
        h3po4_loss=h3po4_loss,
        custom_loss_str=custom_loss_str,
        stateful=stateful
    )

