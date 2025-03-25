import urllib

import requests
import peptacular as pt

import constants

import dataclasses
from typing import Set, Tuple, List, Dict

from msms_compression import SpectrumCompressorF32LzstringUri, SpectrumCompressorUrl

from msms_compression import BaseCompressor, url_encoder, SpectrumCompressorI32, brotli_compressor

lossy_compressor = BaseCompressor(SpectrumCompressorI32(2, 1), brotli_compressor, url_encoder)


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
    neutral_losses: List[str]
    custom_losses: List[Tuple[str, float]]
    isotopes: int
    filter_missing_mono: bool
    filter_interrupted_iso: bool
    min_mz: float
    max_mz: float
    compression_algorithm: str
    min_intensity_type: str
    min_charge: int
    max_charge: int
    top_n: int
    bottom_n: int
    immonium_ions: bool
    internal_fragment_types: Set[str]


class InvalidQueryParam(Exception):
    pass


def serialize_fragments(fragments) -> str:
    return ';'.join([f'{f.count("+")}{f[f.count("+"):]}' for f in fragments])


def deserialize_fragments(s) -> set:
    if not s:
        return set()
    return set('+' * int(f[0]) + f[1:] for f in s.split(';'))


def serialize_spectra(spectra: List[Tuple[float, float]], compression_algorithm: str) -> str:
    if not spectra:
        mzs = []
        ints = []
    else:
        mzs, ints = zip(*spectra)

    if compression_algorithm == 'lzstring':
        return SpectrumCompressorF32LzstringUri.compress(mzs, ints)
    elif compression_algorithm == 'brotli':
        return SpectrumCompressorUrl.compress(mzs, ints)
    elif compression_algorithm == 'lossy':
        return lossy_compressor.compress(mzs, ints)
    elif compression_algorithm == 'key':
        # Assuming 'spectra' is a list with one element that is the key
        data = {"mzs": list(mzs), "intensities": list(ints)}

        response = requests.post(f'{constants.COMP_API}/store', json=data)
        if response.status_code == 200:
            response_data = response.json()
            return response_data.get("key")
        else:
            raise Exception("Failed to store spectra on the server")

    else:
        raise ValueError(f'Invalid compression algorithm: {compression_algorithm}')


def deserialize_spectra(s: str, compression_algorithm: str) -> List[Tuple[float, float]]:
    if compression_algorithm == 'lzstring':
        mzs, ints = SpectrumCompressorF32LzstringUri.decompress(s)
        return list(zip(mzs, ints))
    elif compression_algorithm == 'brotli':
        mzs, ints = SpectrumCompressorUrl.decompress(s)
        return list(zip(mzs, ints))
    elif compression_algorithm == 'lossy':
        mzs, ints = lossy_compressor.decompress(s)
        return list(zip(mzs, ints))
    elif compression_algorithm == 'key':
        # Assuming 'spectra' is a list with one element that is the key
        response = requests.get(f'{constants.COMP_API}/retrieve/{s}')
        if response.status_code == 200:
            response_data = response.json()
            return list(zip(response_data.get("mzs"), response_data.get("intensities")))
        else:
            raise Exception("Failed to retrieve spectra from the server")
    else:
        raise ValueError(f'Invalid compression algorithm: {compression_algorithm}')


def parse_query_params(params) -> QueryParams:
    # Validate Peptide Sequence
    query_peptide_sequence = params.get('sequence', constants.DEFAULT_SEQUENCE)

    try:
        pt.mass(query_peptide_sequence)
    except Exception:
        raise InvalidQueryParam("Invalid peptide sequence.")

    # Mass Type Validation
    query_mass_type = params.get('mass_type', constants.DEFAULT_MASS_TYPE)
    if query_mass_type not in constants.VALID_MASS_TYPES:
        raise InvalidQueryParam("Invalid mass type.")

    # Fragment Types Validation
    query_fragment_types = params.get('fragment_types', serialize_fragments(constants.DEFAULT_FRAGMENT_TYPES))
    try:
        query_fragment_types = deserialize_fragments(query_fragment_types)
    except ValueError:
        raise InvalidQueryParam(f"Invalid fragment types format: '{query_fragment_types}'.")
    if not query_fragment_types.issubset(constants.COLOR_DICT.keys()):
        InvalidQueryParam(f'Invalid fragment types: "{query_fragment_types}".')


    # Fragment Types Validation
    query_internal_fragment_types = params.get('internal_fragment_types', serialize_fragments(constants.DEFAULT_INTERNAL_FRAGMENT_TYPES))
    try:
        query_internal_fragment_types = deserialize_fragments(query_internal_fragment_types)
    except ValueError:
        raise InvalidQueryParam(f'Invalid internal fragment types format: "{query_internal_fragment_types}".')
    if not query_internal_fragment_types.issubset(constants.COLOR_DICT.keys()):
        InvalidQueryParam(f'Invalid internal fragment types: "{query_fragment_types}".')

    # Mass Tolerance Type Validation
    query_mass_tolerance_type = params.get('mass_tolerance_type', constants.DEFAULT_MASS_TOLERANCE_TYPE)
    if query_mass_tolerance_type not in constants.VALID_MASS_TOLERANCE_TYPES:
        InvalidQueryParam("Invalid mass tolerance type.")

    # Mass Tolerance Validation
    query_mass_tolerance = params.get('mass_tolerance', str(constants.DEFAULT_MASS_TOLERANCE))
    try:
        query_mass_tolerance = float(query_mass_tolerance)
        if query_mass_tolerance <= 0:
            InvalidQueryParam("Mass tolerance must be a positive number.")
    except ValueError:
        InvalidQueryParam("Mass tolerance must be a number.")

    # Peak Assignment Validation
    query_peak_assignment = params.get('peak_assignment', constants.DEFAULT_PEAK_ASSIGNMENT)
    if query_peak_assignment not in constants.VALID_PEAK_ASSIGNMENTS:
        InvalidQueryParam("Invalid peak assignment.")

    # Internal Fragments Validation
    query_internal_fragments = params.get('internal_fragments', str(constants.DEFAULT_INTERNAL_FRAGMENTS))
    if query_internal_fragments not in {'True', 'False'}:
        InvalidQueryParam("Internal fragments must be a boolean value.")
    query_internal_fragments = query_internal_fragments == 'True'

    # Minimum Intensity Validation
    query_min_intensity = params.get('min_intensity', str(constants.DEFAULT_MIN_INTENSITY))
    try:
        query_min_intensity = float(query_min_intensity)
        if query_min_intensity < 0:
            InvalidQueryParam("Minimum intensity must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Minimum intensity must be a number.")

    # Y-axis Scale Validation
    query_y_axis_scale = params.get('y_axis_scale', constants.DEFAULT_YAXIS_SCALE)
    if query_y_axis_scale not in constants.VALID_Y_AXIS_SCALES:
        InvalidQueryParam("Invalid y-axis scale.")

    # Hide Unassigned Peaks Validation
    query_hide_unassigned_peaks = params.get('hide_unassigned_peaks', str(constants.DEFAULT_HIDE_UNASSIGNED_PEAKS))
    if query_hide_unassigned_peaks not in {'True', 'False'}:
        InvalidQueryParam("Hide unassigned peaks must be a boolean value.")
    query_hide_unassigned_peaks = query_hide_unassigned_peaks == 'True'

    query_compression_algorithm = params.get('compression_algorithm', constants.DEFAULT_COMPRESSION_ALGORITHM)
    if query_compression_algorithm not in constants.VALID_COMPRESSION_ALGORITHMS:
        InvalidQueryParam(f"Invalid compression algorithm: {query_compression_algorithm}")

    query_spectra = params.get('spectra', serialize_spectra(constants.DEFAULT_SPECTRA, query_compression_algorithm))
    query_spectra = deserialize_spectra(query_spectra, query_compression_algorithm)
    for s in query_spectra:
        if s[0] < 0 or s[1] < 0:
            InvalidQueryParam("Spectra must be non-negative numbers.")

    # Peak Picker Validation
    query_peak_picker = params.get('peak_picker', str(constants.DEFAULT_PEAK_PICKER))
    if query_peak_picker not in {'True', 'False'}:
        InvalidQueryParam("Peak picker must be a boolean value.")
    query_peak_picker = query_peak_picker == 'True'

    # Peak Picker Minimum Intensity Validation
    query_peak_picker_min_intensity = \
        params.get('peak_picker_min_intensity', str(constants.DEFAULT_PEAK_PICKER_MIN_INTENSITY))
    try:
        query_peak_picker_min_intensity = float(query_peak_picker_min_intensity)
        if query_peak_picker_min_intensity < 0:
            InvalidQueryParam("Peak picker minimum intensity must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Peak picker minimum intensity must be a number.")

    # Peak Picker Mass Tolerance Validation
    query_peak_picker_mass_tolerance = \
        params.get('peak_picker_mass_tolerance', str(constants.DEFAULT_PEAK_PICKER_MASS_TOLERANCE))
    try:
        query_peak_picker_mass_tolerance = float(query_peak_picker_mass_tolerance)
        if query_peak_picker_mass_tolerance <= 0:
            InvalidQueryParam("Peak picker mass tolerance must be a positive number.")
    except ValueError:
        InvalidQueryParam("Peak picker mass tolerance must be a number.")

    # Neutral losses validation
    query_neutral_losses = params.get('neutral_losses', ';'.join(constants.DEFAULT_NEUTRAL_LOSSES))
    if query_neutral_losses == '':
        query_neutral_losses = []
    else:
        query_neutral_losses = [nl for nl in query_neutral_losses.split(';') if nl]
    if query_neutral_losses:
        for nl in query_neutral_losses:
            if nl not in constants.NEUTRAL_LOSSES:
                InvalidQueryParam(f"Invalid neutral loss: {nl}")

    query_custom_losses = params.get('custom_losses', ';'.join(constants.DEFAULT_CUSTOM_LOSSES))
    if query_custom_losses == '':
        query_custom_losses = []
    else:
        try:
            query_custom_losses = [(nl.split(':')[0], float(nl.split(':')[1])) for nl in query_custom_losses.split(';') if nl]
        except ValueError:
            InvalidQueryParam("Custom losses must be numbers.")
        if query_custom_losses:
            for nl in query_custom_losses:
                if nl <= 0:
                    InvalidQueryParam(f"Invalid custom loss: {nl}")

    # Isotopes validation
    query_isotopes = params.get('isotopes', str(constants.DEFAULT_ISOTOPES))
    try:
        query_isotopes = int(query_isotopes)
        if query_isotopes < 0:
            InvalidQueryParam("Isotopes must be a non-negative integer.")
    except ValueError:
        InvalidQueryParam("Isotopes must be an integer.")

    # Filter interrupted isotopes validation
    query_filter_missing_mono = params.get('filter_missing_mono', str(constants.DEFAULT_FILTER_MISSING_MONO))
    if query_filter_missing_mono not in {'True', 'False'}:
        InvalidQueryParam("Filter interrupted isotopes must be a boolean value.")
    query_filter_missing_mono = query_filter_missing_mono == 'True'

    # Filter interrupted isotopes validation
    query_filter_interrupted_iso = \
        params.get('filter_interrupted_iso', str(constants.DEFAULT_FILTER_INTERRUPTED_ISO))
    if query_filter_interrupted_iso not in {'True', 'False'}:
        InvalidQueryParam("Filter interrupted isotopes must be a boolean value.")
    query_filter_interrupted_iso = query_filter_interrupted_iso == 'True'

    # validate min/max mz
    query_min_mz = params.get('min_mz', str(constants.DEFAULT_MIN_MZ))
    try:
        query_min_mz = float(query_min_mz)
        if query_min_mz < 0:
            InvalidQueryParam("Minimum m/z must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Minimum m/z must be a number.")

    query_max_mz = params.get('max_mz', str(constants.DEFAULT_MAX_MZ))
    try:
        query_max_mz = float(query_max_mz)
        if query_max_mz < 0:
            InvalidQueryParam("Maximum m/z must be a non-negative number.")
    except ValueError:
        InvalidQueryParam("Maximum m/z must be a number.")


    # validate min/max charge
    query_min_charge = params.get('min_charge', str(constants.DEFAULT_MIN_CHARGE))
    try:
        query_min_charge = int(query_min_charge)
        if query_min_charge < 0:
            InvalidQueryParam("Minimum charge must be a non-negative integer.")
    except ValueError:
        InvalidQueryParam("Minimum charge must be an integer.")

    query_max_charge = params.get('max_charge', str(constants.DEFAULT_MAX_CHARGE))
    try:
        query_max_charge = int(query_max_charge)
        if query_max_charge < 0:
            InvalidQueryParam("Maximum charge must be a non-negative integer.")
    except ValueError:
        InvalidQueryParam("Maximum charge must be an integer.")

    # validate min intensity type
    query_min_intensity_type = params.get('min_intensity_type', str(constants.DEFAULT_MIN_INTENSITY_TYPE))
    if query_min_intensity_type not in {'absolute', 'relative'}:
        InvalidQueryParam("Minimum intensity type must be 'absolute' or 'relative'.")

    # validate top n peaks
    query_top_n_peaks = params.get('top_n_peaks', str(constants.DEFAULT_TOP_N_PEAKS))
    try:
        query_top_n_peaks = int(query_top_n_peaks)
        if query_top_n_peaks < 0:
            InvalidQueryParam("Top N peaks must be a non-negative integer.")

    except ValueError:
        InvalidQueryParam("Top N peaks must be an integer.")

    # validate bottom n peaks
    query_bottom_n_peaks = params.get('bottom_n_peaks', str(constants.DEFAULT_BOTTOM_N_PEAKS))
    try:
        query_bottom_n_peaks = int(query_bottom_n_peaks)
        if query_bottom_n_peaks < 0:
            InvalidQueryParam("Bottom N peaks must be a non-negative integer.")
    except ValueError:
        InvalidQueryParam("Bottom N peaks must be an integer.")

    # validate immonium ions
    query_immonium_ions = params.get('immonium_ions', str(constants.DEFAULT_IMMONIUM_IONS))
    if query_immonium_ions not in {'True', 'False'}:
        InvalidQueryParam("Immonium ions must be a boolean value.")
    query_immonium_ions = query_immonium_ions == 'True'

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
        max_mz=query_max_mz,
        compression_algorithm=query_compression_algorithm,
        min_charge=query_min_charge,
        max_charge=query_max_charge,
        min_intensity_type=query_min_intensity_type,
        top_n=query_top_n_peaks,
        bottom_n=query_bottom_n_peaks,
        immonium_ions=query_immonium_ions,
        internal_fragment_types=query_internal_fragment_types,
    )


def generate_app_url(qp: QueryParams, comp, debug) -> str:
    base_url = constants.BASE_URL

    if debug:
        base_url = 'http://localhost:8501/'

    params = {}
    params['sequence'] = urllib.parse.quote(qp.sequence)
    params['mass_type'] = qp.mass_type
    params['fragment_types'] = serialize_fragments(qp.fragment_types)
    params['mass_tolerance_type'] = str(qp.mass_tolerance_type)
    params['mass_tolerance'] = str(qp.mass_tolerance)
    params['peak_assignment'] = str(qp.peak_assignment)
    params['internal_fragments'] = qp.internal_fragments
    params['min_intensity'] = str(qp.min_intensity)
    params['y_axis_scale'] = str(qp.y_axis_scale)
    params['hide_unassigned_peaks'] = qp.hide_unassigned_peaks
    params['peak_picker'] = qp.peak_picker
    params['peak_picker_min_intensity'] = str(qp.peak_picker_min_intensity)
    params['peak_picker_mass_tolerance'] = str(qp.peak_picker_mass_tolerance)
    params['neutral_losses'] = ';'.join(qp.neutral_losses)
    params['custom_losses'] = ';'.join(f'{nl[0]}:{nl[1]}' for nl in qp.custom_losses)
    params['isotopes'] = qp.isotopes
    params['filter_missing_mono'] = qp.filter_missing_mono
    params['filter_interrupted_iso'] = qp.filter_interrupted_iso
    params['min_mz'] = str(qp.min_mz)
    params['max_mz'] = str(qp.max_mz)
    params['min_charge'] = str(qp.min_charge)
    params['max_charge'] = str(qp.max_charge)
    params['min_intensity_type'] = qp.min_intensity_type
    params['compression_algorithm'] = qp.compression_algorithm
    params['spectra'] = serialize_spectra(qp.spectra, qp.compression_algorithm)
    params['top_n'] = str(qp.top_n)
    params['bottom_n'] = str(qp.bottom_n)
    params['query_immonium_ions'] = qp.immonium_ions
    params['internal_fragment_types'] = serialize_fragments(qp.internal_fragment_types)
    query_string = '&'.join([f'{key}={value}' for key, value in params.items() if value is not None])
    return f'{base_url}?{query_string}'
