#!/usr/bin/env python3
import pandas as pd
import argparse
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def plot_coverage_and_mapping_quality(coverage_df, mapping_quality_df, samplename):
    # Create figure with 2 subplots, top one for coverage, bottom one for mapping quality
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('Genomic Read Depth Coverage', 'Mapping Quality Across Reference'),
        vertical_spacing=0.12,
        row_heights=[0.5, 0.5]
    )

    # Add the first plot for coverage
    fig.add_trace(go.Scatter(
        x=coverage_df['Position (bp)'],
        y=coverage_df['Coverage'],
        mode='lines',
        name='Coverage',
        line=dict(color='#FF6B35', width=2.5),  # Vibrant orange-red
        hovertemplate='<b>Position:</b> %{x} bp<br>' +
                    '<b>Coverage:</b> %{y:.2f}x<br>' +
                    '<extra></extra>',
        legendgroup='coverage'
    ), row=1, col=1)

    # Add shaded region for standard deviation like qualimap tries to do
    fig.add_trace(go.Scatter(
        x=coverage_df['Position (bp)'],
        y=coverage_df['Coverage'] + coverage_df['Std'],
        mode='lines',
        name='Std Upper',
        line=dict(width=0),
        showlegend=False,
        hoverinfo='skip',
        legendgroup='coverage'
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=coverage_df['Position (bp)'],
        y=coverage_df['Coverage'] - coverage_df['Std'],
        mode='lines',
        name='±1 Std Dev',
        line=dict(width=0),
        fillcolor='rgba(255, 107, 53, 0.25)',
        fill='tonexty',
        showlegend=True,
        hoverinfo='skip',
        legendgroup='coverage'
    ), row=1, col=1)

    # Create panel for mapping quality
    fig.add_trace(go.Scatter(
        x=mapping_quality_df['Position (bp)'],
        y=mapping_quality_df['mapping quality'],
        mode='lines',
        name='Mapping Quality',
        line=dict(color='#9B59B6', width=2.5), 
        hovertemplate='<b>Position:</b> %{x} bp<br>' +
                      '<b>Mapping Quality:</b> %{y:.2f}<br>' +
                      '<extra></extra>',
        legendgroup='mq'
    ), row=2, col=1)

    # Add reference line at max mapping quality
    max_mq = mapping_quality_df['mapping quality'].max()
    fig.add_hline(y=max_mq, line_color="#E74C3C", line_width=2,
                  annotation_text=f"Max MQ",
                  annotation_position="right",
                  row=2, col=1)

    # Scale toggle and make pretty
    fig.update_layout(
        hovermode='x unified',
        plot_bgcolor='#F7FAFC',
        height=900,
        template='plotly_white',
        showlegend=True,
        legend=dict(
            orientation='h',
            yanchor='top',
            y=-0.08,
            xanchor='center',
            x=0.5,
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='#CBD5E0',
            borderwidth=1
        ),
        # Add scale toggle buttons for coverage plot
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                buttons=[
                    dict(
                        args=[{
                            "yaxis.type": "linear",
                            "yaxis.tickmode": "auto"
                        }],
                        label="Linear Scale",
                        method="relayout"
                    ),
                    dict(
                        args=[{
                            "yaxis.type": "log",
                            "yaxis.tickmode": "array",
                            "yaxis.tickvals": [0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000],
                            "yaxis.ticktext": ["0.1", "0.5", "1", "5", "10", "50", "100", "500", "1000", "5000", "10000", "50000", "100000"]
                        }],
                        label="Natural Log (ln) Scale",
                        method="relayout"
                    )
                ],
                pad={"r": 10, "t": 10, "b": 10, "l": 10},
                showactive=True,
                x=1.0,
                xanchor="right",
                y=1.02,
                yanchor="bottom",
                bgcolor='#FFE5DB',
                active=0,
                font=dict(color='#2D3748', size=12, family='Arial, sans-serif'),
                bordercolor='#FF6B35',
                borderwidth=2
            )
        ]
    )

    # Update axes for subplot 1 (Coverage)
    fig.update_xaxes(
        title_text='Position (bp)',
        showgrid=True,
        gridcolor='#E2E8F0',
        row=1, col=1
    )
    fig.update_yaxes(
        title_text='Coverage Depth',
        showgrid=True,
        gridcolor='#E2E8F0',
        row=1, col=1
    )

    # Update axes for subplot 2 (Mapping Quality)
    fig.update_xaxes(
        title_text='Position (bp)',
        showgrid=True,
        gridcolor='#E2E8F0',
        row=2, col=1
    )
    fig.update_yaxes(
        title_text='Mapping Quality',
        showgrid=True,
        gridcolor='#E2E8F0',
        range=[0, max_mq * 1.05],
        row=2, col=1
    )

    # Save as HTML with 2 panels
    output_file = f'{samplename}_coverage_plot.html'
    fig.write_html(output_file)
    print(f"Interactive plot saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot genomic read depth coverage and mapping quality.')
    parser.add_argument('coverage_depth_file', help='Path to the coverage data file (TSV format).')
    parser.add_argument('mapping_quality_file', help='Path to the mapping quality data file (TSV format).')
    parser.add_argument('samplename', help='Sample name for output file naming.')
    args = parser.parse_args()

    # Read the coverage data from qualimap output
    coverage_depth_df = pd.read_csv(args.coverage_depth_file, sep='\t')
    coverage_depth_df.columns = coverage_depth_df.columns.str.replace('#', '').str.strip()

    # Read the mapping quality data from qualimap output
    mapping_quality_df = pd.read_csv(args.mapping_quality_file, sep='\t')
    mapping_quality_df.columns = mapping_quality_df.columns.str.replace('#', '').str.strip()

    plot_coverage_and_mapping_quality(coverage_depth_df, mapping_quality_df, args.samplename)