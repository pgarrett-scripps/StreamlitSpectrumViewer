import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def coverage_string(protein_cov_arr, stripped_protein_sequence, c='grey'):
    # Find the maximum coverage value
    max_coverage = max(protein_cov_arr)

    # Color all covered amino acids based on coverage and show index on hover, using a monospace font
    protein_coverage = '<span style="font-family: Courier New, monospace; font-size: 18px;">'
    for i, aa in enumerate(stripped_protein_sequence):
        coverage = protein_cov_arr[i]
        if coverage > 0 and max_coverage > 0:

            protein_coverage += f'<span title="Index: {i + 1}" style="background-color:#e0e0ff; color:{c}; font' \
                                f'-weight:900; padding:3px; margin:1px; border:1px solid #a0a0ff; ' \
                                f'border-radius:3px;">{aa}</span>'
        else:
            protein_coverage += f'<span title="Index: {i + 1}" style="background-color:#f0f0f0; color:#333; ' \
                                f'font-weight:900; padding:3px; margin:1px; border:1px solid #cccccc; ' \
                                f'border-radius:3px;">{aa}</span>'
    protein_coverage += '</span>'

    return protein_coverage

def generate_annonated_spectra_plotly(df, scale='linear', error_scale='ppm'):
    unique_color_labels = df['color_label'].unique().tolist()
    def order(x):
        if x == 'unassigned':
            return 0
        elif x.endswith('i'):
            return (ord('z')+1) * 50 + x.count('+')
        else:
            return ord(x[-1]) * 50 + x.count('+')

    # Move 'unassigned' (or your specific label) to the front
    if 'unassigned' in unique_color_labels:
        unique_color_labels.remove('unassigned')
        unique_color_labels.insert(0, 'unassigned')

    fig_spectra = go.Figure()
    for color_label in unique_color_labels:
        tmp_df = df[df['color_label'] == color_label]
        hover_texts = tmp_df.apply(lambda
                                       row: f"Charge: {row['charge']}<br>M/Z: {row['mz']}<br>Error: {row['error']}"
                                            f"<br>Sequence: {row['sequence']}<br>Label: {row['label']}"
                                            f"<br>Intensity: {row['intensity']}<br>Isotope: {row['isotope']}"
                                            f"<br>Internal: {row['internal']}<br>Ion Type: {row['ion_type']}"
                                            f"<br>loss: {row['loss']}<br>Fragment M/Z: {row['theo_mz']}",
                                   axis=1)

        # Create the Spectra Plot
        fig_spectra.add_trace(go.Bar(x=tmp_df['mz'],
                                     y=tmp_df['intensity'],
                                     width=0.25,
                                     marker_color=tmp_df['color'],
                                     hovertext=hover_texts,
                                     hoverinfo='text',
                                     legendgroup=color_label,
                                     name=color_label,
                                     legendrank=order(color_label),
                                     showlegend=True,
                                     opacity=0.66 if color_label == 'unassigned' else 1.0))


    fig_spectra.update_layout(title='Spectra Plot',
                              xaxis_title='M/Z',
                              yaxis_title='Intensity',
                              yaxis_type=scale)

    fig_error = go.Figure()
    for color_label in df['color_label'].unique():
        if color_label == 'unassigned':
            continue

        tmp_df = df[df['color_label'] == color_label]
        text_colors = tmp_df['color']

        hover_texts = tmp_df.apply(lambda
                                       row: f"Charge: {row['charge']}<br>M/Z: {row['mz']}<br>Error: {row['error']}"
                                            f"<br>Sequence: {row['sequence']}<br>Label: {row['label']}"
                                            f"<br>Intensity: {row['intensity']}<br>Isotope: {row['isotope']}"
                                            f"<br>Internal: {row['internal']}<br>Ion Type: {row['ion_type']}"
                                            f"<br>loss: {row['loss']}<br>Fragment M/Z: {row['theo_mz']}",
                                   axis=1)

        fig_error.add_trace(go.Scatter(x=tmp_df['mz'],
                                       y=tmp_df['error'] if error_scale == 'th' else tmp_df['error_ppm'],
                                       mode='markers+text',
                                       marker_color=tmp_df['color'],
                                       text=tmp_df['label'],
                                       hovertext=hover_texts,
                                       textfont=dict(color=text_colors),
                                       hoverinfo='text',
                                       legendgroup=color_label,
                                       legendrank=order(color_label),
                                       name=color_label,
                                       showlegend=False,
                                       ))

    # Calculate the 95th percentile of the absolute values of the error
    try:
        if error_scale == 'th':
            percentile_95 = np.percentile(df[df['color_label'] != 'unassigned']['error'].abs(), 95)
            mean_error = df[df['color_label'] != 'unassigned']['error'].mean()
            min_error = df[df['color_label'] != 'unassigned']['error'].min()
            max_error = df[df['color_label'] != 'unassigned']['error'].max()
        else:
            percentile_95 = np.percentile(df[df['color_label'] != 'unassigned']['error_ppm'].abs(), 95)
            mean_error = df[df['color_label'] != 'unassigned']['error_ppm'].mean()
            min_error = df[df['color_label'] != 'unassigned']['error_ppm'].min()
            max_error = df[df['color_label'] != 'unassigned']['error_ppm'].max()
    except IndexError:
        percentile_95 = 0
        mean_error = 0
        min_error = -1
        max_error = 1

    # The positive and negative 95th percentiles
    error_positive_95th = mean_error + percentile_95
    error_negative_95th = mean_error - percentile_95

    fig_error.update_layout(title='Error Plot',
                            xaxis_title='M/Z',
                            yaxis_title='Error (th)' if error_scale == 'th' else 'Error (ppm)',
                            yaxis_type=scale)

    fig_error.update_traces(textposition='top center')

    # Combine plots into subplots
    fig = make_subplots(rows=2,
                        cols=1,
                        shared_xaxes=True,
                        vertical_spacing=0.05,
                        row_heights=[1, 2],
                        x_title='M/Z')  # Adjust the relative heights here

    for trace in fig_spectra.data:
        fig.add_trace(trace, row=2, col=1)

    for trace in fig_error.data:
        fig.add_trace(trace, row=1, col=1)

    # Specify the size of the figure
    fig.update_layout(
        width=800,  # Width of the figure in pixels
        height=600,  # Height of the figure in pixels
        title_text="Annotated Spectra"
    )

    # Add horizontal line for positive 95th percentile
    fig.add_shape(
        type='line',
        line=dict(dash='dash', color='grey'),
        x0=df['mz'].min(),
        x1=df['mz'].max(),
        y0=error_positive_95th,
        y1=error_positive_95th,
        row=1,
        col=1,
        opacity=0.5
    )

    # Add horizontal line for negative 95th percentile
    fig.add_shape(
        type='line',
        line=dict(dash='dash', color='grey'),
        x0=df['mz'].min(),
        x1=df['mz'].max(),
        y0=error_negative_95th,
        y1=error_negative_95th,
        row=1,
        col=1,
        opacity=0.5
    )

    fig.add_shape(
        type='line',
        line=dict(dash='dash', color='grey'),
        x0=df['mz'].min(),
        x1=df['mz'].max(),
        y0=mean_error,
        y1=mean_error,
        row=1,
        col=1,
        opacity=0.5
    )

    # Add annotation for positive 95th percentile line
    fig.add_annotation(
        x=df['mz'].max(),
        y=error_positive_95th,
        text=f"+95th percentile ({error_positive_95th:.2f})",
        showarrow=False,
        yshift=10,
        row=1, col=1
    )

    # Add annotation for negative 95th percentile line
    fig.add_annotation(
        x=df['mz'].max(),
        y=error_negative_95th,
        text=f"-95th percentile ({error_negative_95th:.2f})",
        showarrow=False,
        yshift=-10,
        row=1, col=1
    )

    # Add annotation for negative 95th percentile line
    fig.add_annotation(
        x=df['mz'].max(),
        y=mean_error,
        text=f"Mean Error ({mean_error:.2f})",
        showarrow=False,
        yshift=-10,
        row=1, col=1
    )

    # Adding annotations based on ion_type
    for i, row in df.iterrows():
        if row['ion_type']:  # Add annotations only if ion_type is not empty or None
            fig.add_annotation(
                x=row['mz'],
                y=row['intensity'],
                text=row['label'],
                showarrow=False,
                yshift=10,
                font=dict(
                    size=13,
                    color=row['color']
                ),
                row=2, col=1
            )

    # Set y-axis title for the first subplot (row 1)
    fig.update_yaxes(title_text='Mass Error (th)' if error_scale == 'th' else ' Mass Error (ppm)', row=1, col=1)

    # Set y-axis title for the second subplot (row 2)
    fig.update_yaxes(title_text='Intensity', row=2, col=1)
    # Update x-axis for both subplots
    fig.update_yaxes(row=1, col=1, range=[min_error - abs(min_error * 0.3), max_error + abs(max_error * 0.3)])

    return fig


def generate_error_histogram(df, error_scale='ppm'):
    # Select the error column based on the error scale
    error_column = 'error' if error_scale == 'th' else 'error_ppm'

    df = df[df['color_label'] != 'unassigned']

    # Create the histogram
    fig_histogram = go.Figure()
    fig_histogram.add_trace(
        go.Histogram(
            x=df[error_column],
            nbinsx=20,
            marker_color='blue',  # You can choose a different color
        )
    )

    # Update the layout
    fig_histogram.update_layout(
        title='Fragment Ion Error Histogram',
        xaxis_title=f'Error ({error_scale})',
        yaxis_title='Count',
        bargap=0.1,  # Gap between bars
    )

    return fig_histogram


def generate_fragment_plot(unmodified_sequence, spectra_df, internal=False):
    # create a histogram of the internal fragments starting at the N-terminus
    peptide_length = len(unmodified_sequence)
    amino_acid_labels = {i: aa for i, aa in enumerate(unmodified_sequence, start=0)}

    if internal is True:
        spectra_df = spectra_df[spectra_df['internal'] == True]
        reverse_internal_df = spectra_df[spectra_df['ion_type'].str.contains('x|y|z')]
        forward_internal_df = spectra_df[spectra_df['ion_type'].str.contains('a|b|c')]

    else:
        reverse_internal_df = spectra_df[spectra_df['ion_type'].str.contains('a|b|c')]
        forward_internal_df = spectra_df[spectra_df['ion_type'].str.contains('x|y|z')]

    # Aggregate intensities for each amino acid position
    reverse_intensity_agg = reverse_internal_df.groupby('end')['intensity'].sum()
    forward_intensity_agg = forward_internal_df.groupby('start')['intensity'].sum()

    # Create a subplot figure with 2 rows
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Internal Fragment Count", "Internal Fragment Intensity"))

    # Histogram for Fragment Count
    fig.add_trace(go.Histogram(
        x=reverse_internal_df['end'],
        name='Reverse',
        marker_color='red',
        xbins=dict(size=0.5)
    ), row=1, col=1)

    fig.add_trace(go.Histogram(
        x=forward_internal_df['start'],
        name='Forward',
        marker_color='blue',
        xbins=dict(size=0.5)
    ), row=1, col=1)

    # Bar plot for Fragment Intensity
    fig.add_trace(go.Bar(
        x=list(reverse_intensity_agg.index),
        y=reverse_intensity_agg,
        width=0.5,
        name='Reverse Intensity',
        marker_color='red'
    ), row=2, col=1)
    fig.add_trace(go.Bar(
        x=list(forward_intensity_agg.index),
        y=forward_intensity_agg,
        width=0.5,
        name='Forward Intensity',
        marker_color='blue'
    ), row=2, col=1)

    # Update layout and axis titles
    fig.update_layout(title='Internal Fragments Analysis', showlegend=True, barmode='stack')
    fig.update_xaxes(range=[0, peptide_length], tickvals=list(amino_acid_labels.keys()),
                     ticktext=list(amino_acid_labels.values()), row=1, col=1)
    fig.update_xaxes(range=[0, peptide_length], tickvals=list(amino_acid_labels.keys()),
                     ticktext=list(amino_acid_labels.values()), row=2, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_yaxes(title_text="Intensity", row=2, col=1)

    x_tick_vals = np.arange(0.5, peptide_length + 0.5, 1.0)
    fig.update_xaxes(tickvals=x_tick_vals, row=1, col=1, )
    fig.update_xaxes(tickvals=x_tick_vals, row=2, col=1)

    return fig


def generate_fragment_plot_ion_type(unmodified_sequence, spectra_df):
    # create a histogram of the internal fragments starting at the N-terminus
    peptide_length = len(unmodified_sequence)
    amino_acid_labels = {i: aa for i, aa in enumerate(unmodified_sequence, start=0)}

    spectra_df = spectra_df[spectra_df['internal'] == False]
    spectra_df = spectra_df[spectra_df['ion_type'].str.contains('a|b|c|x|y|z')]

    # sort spectra_df by ion_type
    spectra_df = spectra_df.sort_values(by=['ion_type', 'color_label'])

    fig = make_subplots(rows=2, cols=1, subplot_titles=("Fragmentation Site Count", "Fragmentation Site Intensity"))

    unique_legends = set()
    # Add individual traces for each ion type
    for label in spectra_df['label'].unique():
        ion_df = spectra_df[spectra_df['label'] == label]
        color = ion_df['color'].iloc[0]
        color_label = ion_df['color_label'].iloc[0]

        # Bar plot for Fragment Intensity
        reverse_intensity_agg = ion_df.groupby('start')['intensity'].sum()
        forward_intensity_agg = ion_df.groupby('end')['intensity'].sum()

        if 'a' in label or 'b' in label or 'c' in label:
            show_legend = color_label not in unique_legends
            unique_legends.add(color_label)

            fig.add_trace(go.Histogram(
                x=ion_df['end'],
                name=color_label,
                marker_color=color,
                legendgroup=color_label,
                xbins=dict(size=0.5),
                showlegend=show_legend
            ), row=1, col=1)

            fig.add_trace(go.Bar(
                x=list(forward_intensity_agg.index),
                y=forward_intensity_agg,
                name=color_label,
                marker_color=color,
                legendgroup=color_label,
                width=0.5,
                showlegend=False  # Set to False for subsequent traces
            ), row=2, col=1)

        else:
            show_legend = color_label not in unique_legends
            unique_legends.add(color_label)

            fig.add_trace(go.Histogram(
                x=ion_df['start'],
                name=color_label,
                marker_color=color,
                legendgroup=color_label,
                xbins=dict(size=0.5),
                showlegend=show_legend
            ), row=1, col=1)

            fig.add_trace(go.Bar(
                x=list(reverse_intensity_agg.index),
                y=reverse_intensity_agg,
                name=color_label,
                marker_color=color,
                legendgroup=color_label,
                width=0.5,
                showlegend=False  # Set to False for subsequent traces
            ), row=2, col=1)


    # Update layout and axis titles
    fig.update_layout(title='Fragment Analysis', showlegend=True, barmode='stack')
    fig.update_xaxes(range=[0, peptide_length], tickvals=list(amino_acid_labels.keys()),
                     ticktext=list(amino_acid_labels.values()), row=1, col=1)
    fig.update_xaxes(range=[0, peptide_length], tickvals=list(amino_acid_labels.keys()),
                     ticktext=list(amino_acid_labels.values()), row=2, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)
    fig.update_yaxes(title_text="Intensity", row=2, col=1)

    x_tick_vals = np.arange(0.5, peptide_length + 0.5, 1.0)
    fig.update_xaxes(tickvals=x_tick_vals, row=1, col=1)
    fig.update_xaxes(tickvals=x_tick_vals, row=2, col=1)

    # Define the min and max for the x-axis
    x_min = -1  # minimum value for x-axis
    x_max = peptide_length + 1 # maximum value for x-axis, assuming peptide_length is the length of your sequence


    # Update x-axis for both subplots
    fig.update_xaxes(range=[x_min, x_max], row=1, col=1)
    fig.update_xaxes(range=[x_min, x_max], row=2, col=1)

    return fig