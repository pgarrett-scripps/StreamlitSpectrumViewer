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


def generate_annonated_spectra_plotly(df, scale='linear'):
    fig_spectra = go.Figure()
    for ion_type in df['color_label'].unique():
        tmp_df = df[df['color_label'] == ion_type]
        hover_texts = tmp_df.apply(lambda
                                       row: f"Charge: {row['charge']}<br>M/Z: {row['mz']}<br>Error: {row['error']}"
                                            f"<br>Sequence: {row['sequence']}<br>Label: {row['label']}"
                                            f"<br>Intensity: {row['intensity']}<br>Isotope: {row['isotope']}"
                                            f"<br>Internal: {row['internal']}<br>Ion Type: {row['ion_type']}"
                                            f"<br>loss: {row['loss']}",
                                   axis=1)

        # Create the Spectra Plot
        fig_spectra.add_trace(go.Bar(x=tmp_df['mz'],
                                     y=tmp_df['intensity'],
                                     marker_color=tmp_df['color'],
                                     hovertext=hover_texts,
                                     hoverinfo='text',
                                     legendgroup=ion_type,
                                     name=ion_type,
                                     legendrank=ord(ion_type[-1]) * 3 + ion_type.count(
                                         '+') if ion_type != 'unassigned' else 0,
                                     showlegend=True,
                                     opacity=0.7))

    fig_spectra.update_layout(title='Spectra Plot',
                              xaxis_title='M/Z',
                              yaxis_title='Intensity',
                              yaxis_type=scale)

    fig_error = go.Figure()
    for ion_type in df['color_label'].unique():
        if ion_type == 'unassigned':
            continue

        tmp_df = df[df['color_label'] == ion_type]
        text_colors = tmp_df['color']

        hover_texts = tmp_df.apply(lambda
                                       row: f"Charge: {row['charge']}<br>M/Z: {row['mz']}<br>Error: {row['error']}"
                                            f"<br>Sequence: {row['sequence']}<br>Label: {row['label']}"
                                            f"<br>Intensity: {row['intensity']}<br>Isotope: {row['isotope']}"
                                            f"<br>Internal: {row['internal']}<br>Ion Type: {row['ion_type']}"
                                            f"<br>loss: {row['loss']}",
                                   axis=1)

        fig_error.add_trace(go.Scatter(x=tmp_df['mz'],
                                       y=tmp_df['error'],
                                       mode='markers+text',
                                       marker_color=tmp_df['color'],
                                       text=tmp_df['label'],
                                       hovertext=hover_texts,
                                       textfont=dict(color=text_colors),
                                       hoverinfo='text',
                                       legendgroup=ion_type,
                                       legendrank=ord(ion_type[-1]) * 3 + ion_type.count(
                                           '+') if ion_type != 'unassigned' else 0,
                                       name=ion_type,
                                       showlegend=False,
                                       ))

    # Calculate the 95th percentile of the absolute values of the error
    try:
        percentile_95 = np.percentile(df[df['color_label'] != 'unassigned']['error'].abs(), 95)
        mean_error = df[df['color_label'] != 'unassigned']['error'].mean()
    except IndexError:
        percentile_95 = 0
        mean_error = 0

    # The positive and negative 95th percentiles
    error_positive_95th = mean_error + percentile_95
    error_negative_95th = mean_error - percentile_95

    fig_error.update_layout(title='Error Plot',
                            xaxis_title='M/Z',
                            yaxis_title='Error',
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
    fig.update_yaxes(title_text='Mass Error', row=1, col=1)

    # Set y-axis title for the second subplot (row 2)
    fig.update_yaxes(title_text='Intensity', row=2, col=1)

    return fig
