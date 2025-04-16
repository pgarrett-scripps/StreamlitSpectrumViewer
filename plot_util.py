import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import peptacular as pt


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


def generate_annonated_spectra_plotly(df, scale='linear', 
                                      error_scale='ppm', 
                                      line_width=0.25, 
                                      text_size=13, 
                                      marker_size=5,
                                      axis_text_size=12,
                                      title_text_size=20,
                                      tick_text_size=10,
                                      fig_width=1200,
                                      fig_height=800,
                                      hide_error_precentile_labels=False,
                                      hide_error_labels=True,
                                      bold_labels=True):

    df = df.copy(deep=True)

    def format_label(row):
        if row['ion_group_label'] == 'unassigned':
            return ''

        charge_str = "+" + str(row["charge"])
        ion_type_str = row["ion_type"]
        ion_number = row["number"]
        isotope_number = row["isotope"]
        loss = row["loss"]

        ion_label = f"<sup>{charge_str}</sup>{ion_type_str}<sub>{ion_number}</sub>"

        if isotope_number != 0:
            ion_label += f"<sub>[{isotope_number}]</sub>"

        if loss != 0:
            ion_label += f"<sub>({str(round(row['loss'], 2))})</sub>"

        return ion_label
    
    def format_group_label(row):
        if row['ion_group_label'] == 'unassigned':
            return 'unassigned'

        ion_type_str = row["ion_type"]
        charge_str = "+" + str(row["charge"])

        ion_label = f"<sup>{charge_str}</sup>{ion_type_str}"
        return ion_label

    df['format_label'] = df.apply(format_label, axis=1)
    df['format_group_label'] = df.apply(format_group_label, axis=1)

    if bold_labels:
        df['format_label'] = df['format_label'].apply(lambda x: f"<b>{x}</b>")

    unique_color_labels = df['ion_group_label'].unique().tolist()

    def order(x):
        if x == 'unassigned':
            return 0
        elif x.endswith('i'):
            return (ord('z') + 1) * 50 + x.count('+')
        else:
            return ord(x[-1]) * 50 + x.count('+')

    # Move 'unassigned' (or your specific label) to the front
    if 'unassigned' in unique_color_labels:
        unique_color_labels.remove('unassigned')
        unique_color_labels.insert(0, 'unassigned')

    fig_spectra = go.Figure()
    for color_label in unique_color_labels:
        tmp_df = df[df['ion_group_label'] == color_label]
        format_group_label = tmp_df['format_group_label'].iloc[0]
        
        if color_label == 'unassigned':
            hover_texts = tmp_df.apply(lambda
                                        row: f"m/z: {row['mz']}<br>Intensity: {row['intensity']}",
                                    axis=1)
        else:
            hover_texts = tmp_df.apply(lambda
                                        row: f"Charge: {row['charge']}<br>M/Z: {row['mz']}<br>Error: {row['error']}"
                                                f"<br>Sequence: {row['sequence']}<br>Label: {row['label']}"
                                                f"<br>Intensity: {row['intensity']}<br>Isotope: {row['isotope']}"
                                                f"Ion Type: {row['ion_type']}"
                                                f"<br>loss: {row['loss']}<br>Fragment M/Z: {row['theo_mz']}",
                                    axis=1)
        
        # First add all the bar lines
        first = True
        for i, row in tmp_df.iterrows():
            fig_spectra.add_trace(go.Scattergl(
            x=[row['mz'], row['mz']],
            y=[0, row['intensity']],
            mode='lines',
            line=dict(
                width=line_width,
                color=row['color']
            ),
            name=format_group_label,
            legendgroup=color_label,
            legendrank=order(color_label),
            showlegend=first,
            opacity=0.66 if color_label == 'unassigned' else 1.0
            ))
            first = False


        fig_spectra.add_trace(go.Scatter(
            x=tmp_df['mz'],
            y=tmp_df['intensity'],
            mode='markers+text',
            marker=dict(
            size=line_width,  # Invisible markers, just for positioning annotations
            color=tmp_df['color']
            ),
            text=tmp_df['format_label'],
            textposition='top center',
            textfont=dict(
                size=text_size,  # Use textsize parameter
                color=tmp_df['color'],
            ),
            hovertext=hover_texts,
            hoverinfo='text',
            legendgroup=color_label,
            legendrank=order(color_label),
            name=format_group_label,
            showlegend=False,

        ))

    fig_spectra.update_layout(title='Spectra Plot',
                              yaxis_type=scale,
                              title_font=dict(size=title_text_size),
                              xaxis=dict(titlefont=dict(size=axis_text_size),
                                         tickfont=dict(size=axis_text_size/2)),
                              yaxis=dict(titlefont=dict(size=axis_text_size),
                                         tickfont=dict(size=axis_text_size/2)))

    fig_error = go.Figure()
    for color_label in df['ion_group_label'].unique():
        
        if color_label == 'unassigned':
            continue

        tmp_df = df[df['ion_group_label'] == color_label]
        text_colors = tmp_df['color']
        format_group_label = tmp_df['format_group_label'].iloc[0]

        hover_texts = tmp_df.apply(lambda
                                       row: f"Charge: {row['charge']}<br>M/Z: {row['mz']}<br>Error: {row['error']}"
                                            f"<br>Sequence: {row['sequence']}<br>Label: {row['label']}"
                                            f"<br>Intensity: {row['intensity']}<br>Isotope: {row['isotope']}"
                                            f"<br>Ion Type: {row['ion_type']}"
                                            f"<br>loss: {row['loss']}<br>Fragment M/Z: {row['theo_mz']}",
                                   axis=1)

        fig_error.add_trace(go.Scatter(x=tmp_df['mz'],
                                       y=tmp_df['error'] if error_scale == 'th' else tmp_df['error_ppm'],
                                       mode='markers',
                                       marker=dict(
                                           color=tmp_df['color'],
                                           size=marker_size  # Scale marker size based on line_width
                                       ),
                                       text=tmp_df['label'] if not hide_error_labels else '',
                                       hovertext=hover_texts,
                                       textfont=dict(
                                           color=text_colors,
                                           size=text_size  # Use textsize parameter
                                       ),
                                       hoverinfo='text',
                                       legendgroup=color_label,
                                       legendrank=order(color_label),
                                       name=format_group_label,
                                       showlegend=False,
                                       ))

    # Calculate the 95th percentile of the absolute values of the error
    try:
        if error_scale == 'th':
            percentile_95 = np.percentile(df[df['ion_group_label'] != 'unassigned']['error'].abs(), 95)
            mean_error = df[df['ion_group_label'] != 'unassigned']['error'].mean()
            min_error = df[df['ion_group_label'] != 'unassigned']['error'].min()
            max_error = df[df['ion_group_label'] != 'unassigned']['error'].max()
        else:
            percentile_95 = np.percentile(df[df['ion_group_label'] != 'unassigned']['error_ppm'].abs(), 95)
            mean_error = df[df['ion_group_label'] != 'unassigned']['error_ppm'].mean()
            min_error = df[df['ion_group_label'] != 'unassigned']['error_ppm'].min()
            max_error = df[df['ion_group_label'] != 'unassigned']['error_ppm'].max()
    except IndexError:
        percentile_95 = 0
        mean_error = 0
        min_error = -1
        max_error = 1

    # The positive and negative 95th percentiles
    error_positive_95th = mean_error + percentile_95
    error_negative_95th = mean_error - percentile_95
    fig_error.update_layout(title='Error Plot',
                            yaxis_title='Error (th)' if error_scale == 'th' else 'Error (ppm)',
                            title_font=dict(size=title_text_size),
                            xaxis=dict(titlefont=dict(size=axis_text_size),
                                      tickfont=dict(size=axis_text_size/2)),
                            yaxis=dict(titlefont=dict(size=axis_text_size),
                                      tickfont=dict(size=axis_text_size/2)))
    fig_error.update_traces(textposition='top center')

    # Combine plots into subplots
    fig = make_subplots(rows=2,
                        cols=1,
                        shared_xaxes=True,
                        vertical_spacing=0.05,
                        row_heights=[1, 3])  # Adjust the relative heights here

    for trace in fig_spectra.data:
        fig.add_trace(trace, row=2, col=1)

    for trace in fig_error.data:
        fig.add_trace(trace, row=1, col=1)

    # Specify the size of the figure
    fig.update_layout(
        width=fig_width,  # Width of the figure in pixels
        height=fig_height,  # Height of the figure in pixels
        title_text="Annotated Spectra",
        title_font=dict(size=title_text_size),
        font=dict(size=axis_text_size)
    )

    if not hide_error_precentile_labels:

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
            xshift=-50,
            font=dict(size=axis_text_size/2),
            row=1, col=1
        )

        # Add annotation for negative 95th percentile line
        fig.add_annotation(
            x=df['mz'].max(),
            y=error_negative_95th,
            text=f"-95th percentile ({error_negative_95th:.2f})",
            showarrow=False,
            yshift=-10,
            xshift=-50,
            font=dict(size=axis_text_size/2),
            row=1, col=1
        )

        # Add annotation for mean error
        fig.add_annotation(
            x=df['mz'].max(),
            y=mean_error,
            text=f"Mean Error ({mean_error:.2f})",
            showarrow=False,
            yshift=-10,
            xshift=-50,
            font=dict(size=axis_text_size/2),
            row=1, col=1
        )

            
    # Set y-axis title for the first subplot (row 1)
    fig.update_yaxes(title_text='Mass Error (th)' if error_scale == 'th' else ' Mass Error (ppm)', 
                     row=1, col=1, 
                     titlefont=dict(size=axis_text_size),
                     tickfont=dict(size=tick_text_size))

    # Set y-axis title for the second subplot (row 2)
    fig.update_yaxes(title_text='Intensity', 
                     row=2, col=1,
                     titlefont=dict(size=axis_text_size),
                     tickfont=dict(size=tick_text_size))
    
    # update x axis range too
    x_range = df['mz'].max() - df['mz'].min()
    x_range_offset = abs(x_range * 0.05)

    # update x axis - move the title closer to the axis
    fig.update_xaxes(title_text='M/Z',
                     row=2, col=1,
                     titlefont=dict(size=axis_text_size),
                     tickfont=dict(size=tick_text_size),
                    range=[df['mz'].min() - x_range_offset, df['mz'].max() + x_range_offset])  # Reduce standoff to move title closer to axis
                     
    # hide 0 line
    fig.update_yaxes(zeroline=False, row=1, col=1)
    fig.update_yaxes(zeroline=False, row=2, col=1)

    # Place legend in its own dedicated area below the plot
    fig.update_layout(
        legend=dict(
            orientation="h",  # Horizontal orientation
            yanchor="top",
            y=-0.2,  # Place below the plot
            xanchor="center",
            x=0.5,  # Center horizontally
            font=dict(size=tick_text_size),  # Smaller text for the legend
            bgcolor="rgba(255, 255, 255, 0.7)",  # Semi-transparent background
            bordercolor="lightgrey",
            borderwidth=1
        ),
        margin=dict(b=100),  # Add extra margin at the bottom for the legend
        showlegend=True
    )
    
    # the error figure havs the top labels cutoof update top fig layout such that it has a larger y range
    fig.update_yaxes(range=[min_error - (max_error-min_error)*0.2, max_error + (max_error-min_error)*0.4], row=1, col=1)
    fig.update_yaxes(range=[0, df['intensity'].max() * 1.2], row=2, col=1)


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


def get_fragment_match_table_plotly(params, spectra_df, frag_df):
    import pandas as pd
    import plotly.graph_objects as go
    import numpy as np
    
    # Create combined data similar to original function
    combined_data = {"AA": pt.split(params.sequence)}
    
    # Process each ion type and charge
    for ion in params.fragment_types:
        for charge in params.charges:
            ion_df = frag_df.copy()
            ion_df = ion_df[
                (ion_df["ion_type"] == ion)
                & (ion_df["charge"] == charge)
                & (ion_df["internal"] == False)
                & (ion_df["isotope"] == 0)
                & (ion_df["loss"] == 0)
            ]
            ion_df.sort_values(
                by=["start"] if ion in "xyz" else ["end"],
                inplace=True,
                ascending=False if ion in "xyz" else True,
            )

            # Keep only a single number
            ion_df.drop_duplicates(
                subset=["start"] if ion in "xyz" else ["end"], inplace=True
            )

            frags = ion_df["mz"].tolist()

            if ion in "xyz":
                frags = frags[::-1]

            combined_data[f"{charge}{ion}"] = frags

    # Create the combined DataFrame
    combined_df = pd.DataFrame(combined_data)
    
    # Add the numbering columns
    combined_df["# (abc)"] = list(range(1, len(params.unmodified_sequence) + 1))
    combined_df["# (xyz)"] = list(range(1, len(params.unmodified_sequence) + 1))[::-1]

    # Reorder columns just like in the original function
    combined_cols = combined_df.columns.tolist()
    combined_cols.remove("# (abc)")
    combined_cols.remove("# (xyz)")
    combined_cols.remove("AA")
    forward_cols = [col for col in combined_cols if col[-1] in "abc"]
    reverse_cols = [col for col in combined_cols if col[-1] in "xyz"]

    # Sort columns
    forward_cols.sort()
    reverse_cols.sort(reverse=True)

    new_cols = ["# (abc)"] + forward_cols + ["AA"] + reverse_cols + ["# (xyz)"]
    combined_df = combined_df[new_cols]
    
    # Get matched ions for highlighting
    matched_ions = spectra_df[spectra_df["ion_type"] != ""]
    accepted_normal_ions = matched_ions[matched_ions["internal"] == False]["ion_label"].tolist()
    accepted_internal_ions = matched_ions[matched_ions["internal"] == True]["ion_label"].tolist()
    accepted_internal_ions = [ion[:-1] for ion in accepted_internal_ions]
    
    # Format values with 4 decimal places precision
    formatted_df = combined_df.copy()
    for col in combined_df.columns:
        if col not in ["AA", "# (abc)", "# (xyz)"]:
            formatted_df[col] = formatted_df[col].apply(lambda x: f"{x:.4f}" if not pd.isna(x) else "")
    
    # Create arrays for styling each cell
    cells_bg_color = np.full((len(combined_df), len(new_cols)), 'white')
    cells_text_color = np.full((len(combined_df), len(new_cols)), 'black')
    
    # Set header/index column styling
    for row_idx in range(len(combined_df)):
        for col_idx, col_name in enumerate(new_cols):
            if col_name in ["AA", "# (abc)", "# (xyz)"]:
                cells_bg_color[row_idx, col_idx] = 'gainsboro'
    
    # Apply highlighting based on the logic from the original function
    for row_idx in range(len(combined_df)):
        for col_idx, col_name in enumerate(new_cols):
            if col_name not in ["AA", "# (abc)", "# (xyz)"]:
                ion = col_name[-1]
                charge = int(col_name[:-1])
                
                if ion in "abc":
                    ion_number = row_idx + 1
                else:
                    ion_number = len(params.unmodified_sequence) - row_idx
                    
                ion_key = col_name + str(ion_number)
                mz = combined_df.iloc[row_idx][col_name]
                
                if pd.isna(mz) or mz <= params.min_mz or mz >= params.max_mz:
                    cells_bg_color[row_idx][col_idx] = '#BEBEBE'
                    cells_text_color[row_idx][col_idx] = 'black'
                else:
                    if ion_key in accepted_normal_ions:
                        cells_bg_color[row_idx][col_idx] = params.get_color(ion, charge)
                        cells_text_color[row_idx][col_idx] = 'white'
                    elif ion_key in accepted_internal_ions:
                        cells_bg_color[row_idx][col_idx] = params.get_color(ion, charge)
                        cells_text_color[row_idx][col_idx] = 'magenta'
    
    # Create the Plotly table using new approach with cell data
    cells_values = [formatted_df[col].tolist() for col in formatted_df.columns]
    
    # Create figure
    fig = go.Figure()
    
    # Add the table
    fig.add_trace(
        go.Table(
            header=dict(
                values=list(formatted_df.columns),
                fill_color='gainsboro',
                font=dict(color='black', size=12, family="Arial"),
                align='center',
                height=25
            ),
            cells=dict(
                values=cells_values,
                fill_color=cells_bg_color.tolist(),
                font=dict(
                    color=cells_text_color.tolist(),
                    size=11,
                    family="Arial"
                ),
                align='center',
                height=20
            ),
            columnwidth=[30] * len(new_cols)  # Make columns narrower
        )
    )
    
    # Update layout for smaller size and less padding
    fig.update_layout(
        margin=dict(l=5, r=5, t=5, b=5),
        height=25 * (len(formatted_df) + 1),
        width=min(800, 40 * len(formatted_df.columns)),
    )
    
    return fig

# Example usage
# table = get_fragment_match_table_plotly(params, spectra_df, frag_df)
# table.show()


def generate_fragment_plot_ion_type(unmodified_sequence, spectra_df):
    # create a histogram of the internal fragments starting at the N-terminus
    peptide_length = len(unmodified_sequence)
    amino_acid_labels = {i: aa for i, aa in enumerate(unmodified_sequence, start=0)}

    spectra_df = spectra_df[spectra_df['internal'] == False]
    spectra_df = spectra_df[spectra_df['ion_type'].str.contains('a|b|c|x|y|z')]

    # sort spectra_df by ion_type
    spectra_df = spectra_df.sort_values(by=['ion_type', 'ion_group_label'])

    fig = make_subplots(rows=2, cols=1, subplot_titles=("Fragmentation Site Count", "Fragmentation Site Intensity"))

    unique_legends = set()
    # Add individual traces for each ion type
    for label in spectra_df['label'].unique():
        ion_df = spectra_df[spectra_df['label'] == label]
        color = ion_df['color'].iloc[0]
        color_label = ion_df['ion_group_label'].iloc[0]

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
    x_max = peptide_length + 1  # maximum value for x-axis, assuming peptide_length is the length of your sequence

    # Update x-axis for both subplots
    fig.update_xaxes(range=[x_min, x_max], row=1, col=1)
    fig.update_xaxes(range=[x_min, x_max], row=2, col=1)

    return fig
