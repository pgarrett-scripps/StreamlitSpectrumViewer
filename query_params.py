import urllib

from peptacular.mass import valid_mass_sequence

import constants

import dataclasses
from typing import Set, Tuple, List

from lzstring import LZString

def decompress_from_encoded_uri_component(compressed_string):
    decompressed = LZString.decompressFromEncodedURIComponent(compressed_string)
    return decompressed

def compress_uri_component(string):
    compressed = LZString.compressToEncodedURIComponent(string)
    return compressed


@dataclasses.dataclass
class QueryParams:
    sequence: str
    mass_type: str
    fragment_types: Set[str]
    mass_tolerance_type: str
    mass_tolerance: float
    peak_assignment: str
    internal_fragments: bool
    min_intensity: float
    y_axis_scale: str
    hide_unassigned_peaks: bool
    spectra: List[Tuple[float, float]]
    peak_picker: bool
    peak_picker_min_intensity: float
    peak_picker_mass_tolerance: float
    neutral_losses: Set[str]
    custom_losses: Set[float]
    isotopes: int
    filter_missing_mono: bool
    filter_interrupted_iso: bool
    min_mz: float
    max_mz: float

    # Add any other fields that you need

    # Add any additional methods or validations as needed


class InvalidQueryParam(Exception):
    pass


def serialize_fragments(fragments) -> str:
    return ';'.join([f'{f.count("+")}{f[-1]}' for f in fragments])


def deserialize_fragments(s) -> set:
    return set('+' * int(f[0]) + f[1] for f in s.split(';'))


def serialize_spectra(spectra: List[Tuple[float, float]]) -> str:
    return ';'.join([f'{s[0]}:{s[1]}' for s in spectra])


def deserialize_spectra(s: str) -> List[Tuple[float, float]]:
    return [(float(f.split(':')[0]), float(f.split(':')[1])) for f in s.split(';')]


def parse_query_params(params) -> QueryParams:
    # Validate Peptide Sequence
    query_peptide_sequence = params.get('sequence', [constants.DEFAULT_SEQUENCE])[0]
    if valid_mass_sequence(query_peptide_sequence) is False:
        raise InvalidQueryParam("Invalid peptide sequence.")

    # Mass Type Validation
    query_mass_type = params.get('mass_type', [constants.DEFAULT_MASS_TYPE])[0]
    if query_mass_type not in constants.VALID_MASS_TYPES:
        raise InvalidQueryParam("Invalid mass type.")

    # Fragment Types Validation
    query_fragment_types = params.get('fragment_types', [serialize_fragments(constants.DEFAULT_FRAGMENT_TYPES)])[0]
    try:
        query_fragment_types = deserialize_fragments(query_fragment_types)
    except ValueError:
        raise InvalidQueryParam("Invalid fragment types format.")
    if not query_fragment_types.issubset(constants.COLOR_DICT.keys()):
        InvalidQueryParam(f"Invalid fragment types: {query_fragment_types}")

    # Mass Tolerance Type Validation
    query_mass_tolerance_type = params.get('mass_tolerance_type', [constants.DEFAULT_MASS_TOLERANCE_TYPE])[0]
    if query_mass_tolerance_type not in constants.VALID_MASS_TOLERANCE_TYPES:
        InvalidQueryParam("Invalid mass tolerance type.")

    # Mass Tolerance Validation
    query_mass_tolerance = params.get('mass_tolerance', [str(constants.DEFAULT_MASS_TOLERANCE)])[0]
    try:
        query_mass_tolerance = float(query_mass_tolerance)
        if query_mass_tolerance <= 0:
            InvalidQueryParam("Mass tolerance must be a positive number.")
    except ValueError:
        InvalidQueryParam("Mass tolerance must be a number.")

    # Peak Assignment Validation
    query_peak_assignment = params.get('peak_assignment', [constants.DEFAULT_PEAK_ASSIGNMENT])[0]
    if query_peak_assignment not in constants.VALID_PEAK_ASSIGNMENTS:
        InvalidQueryParam("Invalid peak assignment.")

    # Internal Fragments Validation
    query_internal_fragments = params.get('internal_fragments', [str(constants.DEFAULT_INTERNAL_FRAGMENTS)])[0]
    if query_internal_fragments not in {'True', 'False'}:
        InvalidQueryParam("Internal fragments must be a boolean value.")
    query_internal_fragments = query_internal_fragments == 'True'

    # Minimum Intensity Validation
    query_min_intensity = params.get('min_intensity', [str(constants.DEFAULT_MIN_INTENSITY)])[0]
    try:
        query_min_intensity = float(query_min_intensity)
        if query_min_intensity < 0:
            InvalidQueryParam("Minimum intensity must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Minimum intensity must be a number.")

    # Y-axis Scale Validation
    query_y_axis_scale = params.get('y_axis_scale', [constants.DEFAULT_YAXIS_SCALE])[0]
    if query_y_axis_scale not in constants.VALID_Y_AXIS_SCALES:
        InvalidQueryParam("Invalid y-axis scale.")

    # Hide Unassigned Peaks Validation
    query_hide_unassigned_peaks = params.get('hide_unassigned_peaks', [str(constants.DEFAULT_HIDE_UNNASSIGNED_PEAKS)])[
        0]
    if query_hide_unassigned_peaks not in {'True', 'False'}:
        InvalidQueryParam("Hide unassigned peaks must be a boolean value.")
    query_hide_unassigned_peaks = query_hide_unassigned_peaks == 'True'

    query_spectra = params.get('spectra', [compress_uri_component(serialize_spectra(constants.DEFAULT_SPECTRA))])[0]
    query_spectra = deserialize_spectra(decompress_from_encoded_uri_component(query_spectra))
    for s in query_spectra:
        if s[0] < 0 or s[1] < 0:
            InvalidQueryParam("Spectra must be non-negative numbers.")

    # Peak Picker Validation
    query_peak_picker = params.get('peak_picker', [str(constants.DEFAULT_PEAK_PICKER)])[0]
    if query_peak_picker not in {'True', 'False'}:
        InvalidQueryParam("Peak picker must be a boolean value.")
    query_peak_picker = query_peak_picker == 'True'

    # Peak Picker Minimum Intensity Validation
    query_peak_picker_min_intensity = \
        params.get('peak_picker_min_intensity', [str(constants.DEFAULT_PEAK_PICKER_MIN_INTENSITY)])[0]
    try:
        query_peak_picker_min_intensity = float(query_peak_picker_min_intensity)
        if query_peak_picker_min_intensity < 0:
            InvalidQueryParam("Peak picker minimum intensity must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Peak picker minimum intensity must be a number.")

    # Peak Picker Mass Tolerance Validation
    query_peak_picker_mass_tolerance = \
        params.get('peak_picker_mass_tolerance', [str(constants.DEFAULT_PICK_PEAKER_MASS_TOLERANCE)])[0]
    try:
        query_peak_picker_mass_tolerance = float(query_peak_picker_mass_tolerance)
        if query_peak_picker_mass_tolerance <= 0:
            InvalidQueryParam("Peak picker mass tolerance must be a positive number.")
    except ValueError:
        InvalidQueryParam("Peak picker mass tolerance must be a number.")

    # Neutral losses validation
    query_neutral_losses = params.get('neutral_losses', [''])[0]
    query_neutral_losses = set(nl for nl in query_neutral_losses.split(';') if nl)
    if query_neutral_losses:
        for nl in query_neutral_losses:
            if nl not in constants.NEUTRAL_LOSSES:
                InvalidQueryParam(f"Invalid neutral loss: {nl}")

    query_custom_losses = params.get('custom_losses', [''])[0]
    try:
        query_custom_losses = set(float(nl) for nl in query_custom_losses.split(';') if nl)
    except ValueError:
        InvalidQueryParam("Custom losses must be numbers.")
    if query_custom_losses:
        for nl in query_custom_losses:
            if nl <= 0:
                InvalidQueryParam(f"Invalid custom loss: {nl}")


    # Isotopes validation
    query_isotopes = params.get('isotopes', [str(constants.DEFAULT_ISOTOPES)])[0]
    try:
        query_isotopes = int(query_isotopes)
        if query_isotopes < 0:
            InvalidQueryParam("Isotopes must be a non-negative integer.")
    except ValueError:
        InvalidQueryParam("Isotopes must be an integer.")

    # Filter interrupted isotopes validation
    query_filter_missing_mono = params.get('filter_missing_mono', [str(constants.DEFAULT_FILTER_MISSING_MONO)])[0]
    if query_filter_missing_mono not in {'True', 'False'}:
        InvalidQueryParam("Filter interrupted isotopes must be a boolean value.")
    query_filter_missing_mono = query_filter_missing_mono == 'True'

    # Filter interrupted isotopes validation
    query_filter_interrupted_iso = \
        params.get('filter_interrupted_iso', [str(constants.DEFAULT_FILTER_INTERRUPTED_ISO)])[0]
    if query_filter_interrupted_iso not in {'True', 'False'}:
        InvalidQueryParam("Filter interrupted isotopes must be a boolean value.")
    query_filter_interrupted_iso = query_filter_interrupted_iso == 'True'

    # validate min/max mz
    query_min_mz = params.get('min_m/z', [str(constants.DEFAULT_MIN_MZ)])[0]
    try:
        query_min_mz = float(query_min_mz)
        if query_min_mz < 0:
            InvalidQueryParam("Minimum m/z must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Minimum m/z must be a number.")

    query_max_mz = params.get('max_m/z', [str(constants.DEFAULT_MAX_MZ)])[0]
    try:
        query_max_mz = float(query_max_mz)
        if query_max_mz < 0:
            InvalidQueryParam("Maximum m/z must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Maximum m/z must be a number.")

    return QueryParams(
        sequence=query_peptide_sequence,
        mass_type=query_mass_type,
        fragment_types=query_fragment_types,
        mass_tolerance_type=query_mass_tolerance_type,
        mass_tolerance=query_mass_tolerance,
        peak_assignment=query_peak_assignment,
        internal_fragments=query_internal_fragments,
        min_intensity=query_min_intensity,
        y_axis_scale=query_y_axis_scale,
        hide_unassigned_peaks=query_hide_unassigned_peaks,
        spectra=query_spectra,
        peak_picker=query_peak_picker,
        peak_picker_min_intensity=query_peak_picker_min_intensity,
        peak_picker_mass_tolerance=query_peak_picker_mass_tolerance,
        neutral_losses=query_neutral_losses,
        custom_losses=query_custom_losses,
        isotopes=query_isotopes,
        filter_missing_mono=query_filter_missing_mono,
        filter_interrupted_iso=query_filter_interrupted_iso,
        min_mz=query_min_mz,
        max_mz=query_max_mz
    )


def generate_app_url(qp: QueryParams) -> str:
    base_url = constants.BASE_URL
    params = {
        'sequence': urllib.parse.quote(qp.sequence),
        'mass_type': qp.mass_type,
        'fragment_types': serialize_fragments(qp.fragment_types),
        'mass_tolerance_type': str(qp.mass_tolerance_type),
        'mass_tolerance': str(qp.mass_tolerance),
        'peak_assignment': str(qp.peak_assignment),
        'internal_fragments': qp.internal_fragments,
        'min_intensity': str(qp.min_intensity),
        'y_axis_scale': str(qp.y_axis_scale),
        'neutral_losses': ';'.join(qp.neutral_losses),
        'custom_losses': ';'.join(str(nl) for nl in qp.custom_losses),
        'hide_unassigned_peaks': qp.hide_unassigned_peaks,
        'peak_picker': qp.peak_picker,
        'peak_picker_min_intensity': str(qp.peak_picker_min_intensity),
        'peak_picker_mass_tolerance': str(qp.peak_picker_mass_tolerance),
        'isotopes': str(qp.isotopes),
        'filter_missing_mono': str(qp.filter_missing_mono),
        'filter_interrupted_iso': str(qp.filter_interrupted_iso),
        'min_mz': str(qp.min_mz),
        'max_mz': str(qp.max_mz),
        'spectra': compress_uri_component(serialize_spectra(qp.spectra)),
    }
    query_string = '&'.join([f'{key}={value}' for key, value in params.items() if value is not None])
    return f'{base_url}?{query_string}'
