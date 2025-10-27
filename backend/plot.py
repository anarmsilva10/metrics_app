import plotly.graph_objects as go
import argparse
from pathlib import Path
from backend import depth
from log import shared_log

def coverage_plot(cram_path, bed_df, method='run_pysam_coverage', threshold=20, padding=0):
    """
    Generate a coverage plot for the given CRAM/BAM file and BED regions.

    Args:
        cram_path (str): Path to the CRAM/BAM input file.
        bed_df (pd.DataFrame): Filtered BED DataFrame containing regions of interest.
        method (str, optional): Method to extract depth ('run_samtools_depth', 'run_pysam_depth', 'run_pysam_coverage'). Default is 'run_pysam_coverage'.
        threshold (int, optional): Coverage threshold to highlight low-coverage regions. Default is 20.
        padding (int, optional): Padding applied to exon regions in the plot. Default is 0.

    Returns:
        plotly.graph_objects.Figure | None: Plotly figure object with coverage plot.
        Returns None if no depth data or BED regions are available.
    """

    if bed_df.empty:
        shared_log.logger.warning('No regions found.')
        return None

    call =  getattr(depth, method)
    depth_df = call(cram_path, bed_df)

    if depth_df.empty:
        shared_log.logger.warning('No depth data found.')
        return None
   
    # Annotate depth_df with gene/exon info
    depth_df['Gene'] = 'Unknown'
    depth_df['Exon'] = 'Unknown'
    for _, row in bed_df.iterrows():
        mask = (depth_df['POSITION'] >= row['START']) & (depth_df['POSITION'] <= row['END'])
        depth_df.loc[mask, 'Gene'] = row['GENE']
        depth_df.loc[mask, 'Exon'] = row['EXON']

    # Prepare plot data with gaps for non-contiguous positions
    plot_x, plot_y, hover_text = [], [], []
    positions, depths = depth_df['POSITION'].values, depth_df['DEPTH'].values
    gene_info, exon_info = depth_df['Gene'].values, depth_df['Exon'].values

    for i in range(len(positions)):
        plot_x.append(positions[i])
        plot_y.append(depths[i])
        hover_text.append(f"Gene: {gene_info[i]}<br>Exon: {exon_info[i]}<br>Pos: {positions[i]}<br>Depth: {depths[i]}")
        if i < len(positions) - 1 and positions[i + 1] - positions[i] > 1:
            plot_x.append(None)
            plot_y.append(None)
            hover_text.append(None)

    # Main coverage trace
    coverage_trace = go.Scatter(
        x=plot_x,
        y=plot_y,
        mode='lines',
        name='Depth',
        text=hover_text,
        hoverinfo='text',
        line=dict(color='blue'),
        connectgaps=False
    )

    # Mean and threshold lines
    mean_coverage = depth_df['DEPTH'].mean()
    mean_line = go.Scatter(
        x=[min(positions), max(positions)],
        y=[mean_coverage, mean_coverage],
        mode='lines',
        name=f'Average Depth: {mean_coverage:.1f}X',
        line=dict(color='green'),
        connectgaps=False,
        hoverinfo='skip'
    )

    threshold_line = go.Scatter(
        x=[min(positions), max(positions)],
        y=[threshold, threshold],
        mode='lines',
        name=f'Threshold: {threshold}X',
        line=dict(color='red'),
        connectgaps=False,
        hoverinfo='skip'
    )

    # Highlight low-coverage regions
    highlight_traces = []
    highlight_x, highlight_y = [], []
    show_legend_once = True

    for i in range(len(plot_x)):
        x_val, y_val = plot_x[i], plot_y[i]
        if y_val is not None and y_val <= threshold:
            highlight_x.append(x_val)
            highlight_y.append(y_val)
        else:
            if highlight_x:
                highlight_traces.append(
                    go.Scatter(
                        x=highlight_x,
                        y=highlight_y,
                        fill='tozeroy',
                        mode='lines',
                        fillcolor='rgba(255, 0, 0, 0.2)',
                        line=dict(color='rgba(255,0,0,0)'),
                        name='Below Threshold' if show_legend_once else None,
                        showlegend=show_legend_once,
                        connectgaps=False,
                        hoverinfo='skip'
                    )
                )
                highlight_x, highlight_y = [], []
                show_legend_once = False

    if highlight_x:
        highlight_traces.append(
            go.Scatter(
                x=highlight_x,
                y=highlight_y,
                fill='tozeroy',
                mode='lines',
                fillcolor='rgba(255, 0, 0, 0.2)',
                line=dict(color='rgba(255,0,0,0)'),
                name='Below Threshold' if show_legend_once else None,
                showlegend=show_legend_once,
                connectgaps=False,
                hoverinfo='skip'
            )
        )

    # Exon region shapes
    exon_shapes = [
        dict(
            type="rect",
            xref="x",
            yref="paper",
            x0=row['START'],
            x1=row['END'],
            y0=0,
            y1=1,
            fillcolor="LightSkyBlue",
            opacity=0.3,
            layer="below",
            line_width=0
        ) for _, row in bed_df.iterrows()
    ]

    # Fake trace for exon legend
    exon_legend_trace = go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='LightSkyBlue'),
        fillcolor='LightSkyBlue',
        name='Exon',
        showlegend=True,
        hoverinfo='skip'
    )

    # Fake trace for padding legend
    padding_trace = go.Scatter(
        x=[None], y=[None],
        mode='lines',
        line=dict(color='white'),
        name=f'Padding: {padding}',
        showlegend=True,
        hoverinfo='skip'
    )

    cram_name = Path(cram_path).name

    # Assemble all traces and layout
    layout = go.Layout(
        title=f'Depth of Coverage Plot for {cram_name}',
        xaxis=dict(title='', rangeslider=dict(visible=True), tickformat='d', tickangle=45),
        yaxis=dict(title='Depth'),
        hovermode='closest',
        shapes=exon_shapes
    )

    fig = go.Figure(data=[coverage_trace, mean_line, threshold_line, exon_legend_trace, padding_trace] + highlight_traces,
                    layout=layout)

    return fig

def main():
    parser = argparse.ArgumentParser(description="Calculate depth of coverage from CRAM/BAM and BED files.")

    # Mandatory arguments
    parser.add_argument('cram_path', type=str, help="Path to the input CRAM/BAM file.")
    parser.add_argument('bed_path', type=str, help="Path to the input BED file.")

    #  Optional arguments
    parser.add_argument('--log_file', type=str, help='Path to a log file (default: log to stdout).')
    parser.add_argument('--gene_selection', type=str, nargs='+', help="List of genes to include (optional).")
    parser.add_argument('--exon_selection', type=int, nargs='+', help="List of exons to include (optional).")

    # Method argument - It isn't mandatory, it has a default method
    parser.add_argument('--method', choices=['run_samtools_depth', 'run_pysam_depth','run_pysam_coverage'], default='run_pysam_coverage', help="Choose the method to extract depth (default: run_pysam_coverage).")
    parser.add_argument('--output', type=str, help="Path to save the output plot (e.g., plot.html or plot.png)")
    parser.add_argument('--output_format', choices=['html', 'png', 'pdf', 'svg'], default='html', help="Output format (default: html). Use 'png' or 'pdf' for images.")
    parser.add_argument('--threshold', type=int, default=400, help='COverage threshold (default=400).')


    # Parsing arguments
    args = parser.parse_args()

    if args.log_file:
        shared_log.logger = shared_log.Logging(filename=args.log_file, foldername=__path__).get_logger()
   
    shared_log.logger.info(f"Using BED file: {args.bed_path} and BAM/CRAM file: {args.cram_path}.")

    bed_df = depth.filter_bed(args.bed_path, args.gene_selection, args.exon_selection)

    fig = coverage_plot(args.cram_path, bed_df, method=args.method, threshold=args.threshold)

        # Save figure if requested
    if args.output:
        try:
            if args.output_format == 'html':
                fig.write_html(args.output)
            else:
                fig.write_image(args.output)
            shared_log.logger.info(f"Plot saved to {args.output}")
        except Exception as e:
            shared_log.logger.error(f"Failed to save plot: {e}")

# Calling main function
if __name__ == "__main__":
    main()
