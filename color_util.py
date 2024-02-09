import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import colorsys


def convert_color(color, format_type='hex'):
    """
    Convert matplotlib color to various formats.

    :param color: Color in RGBA format
    :param format_type: Format to convert to ('hex', 'rgb', 'hsl', 'hsv')
    :return: Color in the specified format
    """
    if format_type == 'hex':
        return mcolors.to_hex(color)
    elif format_type == 'rgb':
        return f'rgb({int(color[0] * 255)}, {int(color[1] * 255)}, {int(color[2] * 255)})'
    elif format_type == 'hsl':
        h, l, s = colorsys.rgb_to_hls(color[0], color[1], color[2])
        return f'hsl({int(h * 360)}, {int(s * 100)}%, {int(l * 100)}%)'
    elif format_type == 'hsv':
        h, s, v = colorsys.rgb_to_hsv(color[0], color[1], color[2])
        return f'hsv({int(h * 360)}, {int(s * 100)}%, {int(v * 100)}%)'
    else:
        return color  # Default is to return the original color


def get_color_for_state(ion, charge, min_charge, max_charge, format_type='hex'):
    """
    Get color for a given category and state.

    :param ion: Category of the item (e.g., 'i', 'a', etc.)
    :param charge: Current state (number of '+' signs)
    :param max_charge: Maximum state in the dataset
    :return: Color for the given state
    """
    # Define colormaps for each category
    colormaps = {
        'i': plt.cm.spring,
        'a': plt.cm.BrBG,
        'b': plt.cm.winter,
        'c': plt.cm.PRGn,
        'x': plt.cm.PuOr,
        'y': plt.cm.seismic,
        'z': plt.cm.PRGn,
    }

    def get_pos(c):
        return (c - min_charge) / (max_charge - min_charge)

    color_map_lamdas = {
        'i': lambda x: get_pos(x) * 0.4,
        'a': lambda x: (1 - get_pos(x)) * 0.27,
        'b': lambda x: (1 - get_pos(x)) * 0.55,
        'c': lambda x: 0.75 + (get_pos(x)) * 0.25,
        'x': lambda x: (1 - get_pos(x)) * 0.2 + 0.1,
        'y': lambda x: 0.7 + (get_pos(x)) * 0.2,
        'z': lambda x: (1 - get_pos(x)) * 0.3,
    }

    # Get the colormap for the category
    cmap = colormaps.get(ion)  # Default to gray if category not found

    # Calculate the color position
    color_pos = color_map_lamdas.get(ion)(charge)

    # Get the color from the colormap
    return convert_color(cmap(color_pos), format_type=format_type)


def get_color_dict(min_charge, max_charge):
    """
    Get a dictionary of colors for each state.

    :param max_charge: Maximum state in the dataset
    :param min_charge: Minimum state in the dataset
    :return: Dictionary of colors for each state
    """
    ion_to_color = {'unassigned': 'grey'}
    for ion in 'abcxyziI':
        for charge in range(min_charge, max_charge + 1):
            if ion == 'I':
                ion_to_color[f'{"+" * charge}{ion}'] = 'pink'
            else:
                ion_to_color[f'{"+" * charge}{ion}'] = get_color_for_state(ion, charge, min_charge, max_charge)

    return ion_to_color
