import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score 
from scipy.stats import linregress
import geopandas as gpd
from shapely.geometry import Polygon, Point

def get_transcripts_cor(xenium_bundle, sudo_bulk_PATH):
    transcripts_df = pd.read_parquet(str(xenium_bundle + '/transcripts.parquet'))
    sudo_bulk = pd.read_csv(sudo_bulk_PATH)
    
    counts = transcripts_df['feature_name'].value_counts()
    counts_df = counts.reset_index()
    counts_df.columns = ['gene', 'count']

    merged_df = counts_df.merge(sudo_bulk, how='left', on='gene')
    filtered_df = merged_df[
        ~merged_df['gene'].str.startswith(('UnassignedCodeword', 'NegControlCodeword', 'DeprecatedCodeword', 'NegControlProbe'))
    ]
    
    filtered_df = filtered_df[(filtered_df['count'] > 0) & (filtered_df['counts'] > 0)]
    corr_matrix = np.corrcoef(filtered_df['count'], filtered_df['counts'])
    corr = corr_matrix[0,1]
    R_sq = corr**2
    
    x = filtered_df['count']
    y = filtered_df['counts']

    log_x = np.log10(x)
    log_y = np.log10(y)
    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)

    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.7)
    plt.text(
        0.05, 0.95, f"R = {corr:.3f}", 
        transform=plt.gca().transAxes, 
        fontsize=12, 
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)
    )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Xenium')
    plt.ylabel('Sudo-Bulk RNA')
    plt.title('Transcript Correlation')
    #plt.title(sample_name)
    #plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.show()
    
    
def categorize_prefix(value):
    if value.startswith('UnassignedCodeword'):
        return 'UnassignedCodeword'
    elif value.startswith('NegControlProbe'):
        return 'NegControlProbe'
    elif value.startswith('NegControlCodeword'):
        return 'NegControlProbe'
    elif value.startswith('DeprecatedCodeword'):
        return 'DeprecatedCodeword'
    else:
        return 'Probe'

    
def get_transcripts_in_polygon(xenium_bundle, polygon_file):
    
    transcripts_df = pd.read_parquet(str(xenium_bundle + '/transcripts.parquet'))
    transcripts_gdf = gpd.GeoDataFrame(
    transcripts_df,
    geometry=[Point(x, y) for x, y in zip(transcripts_df['x_location'], transcripts_df['y_location'])]
    )
    
    polygon_df = pd.read_csv(polygon_file)
    polygons = (
        polygon_df.groupby('polygon_id')
        .apply(lambda group: Polygon(zip(group['X'], group['Y'])))
        .reset_index(name='geometry')
    )
    polygons_gdf = gpd.GeoDataFrame(polygons, geometry='geometry')
    
    transcripts_in_polygons = gpd.sjoin(transcripts_gdf, polygons_gdf, how='inner', predicate='within')
    
    counts = transcripts_in_polygons['feature_name'].value_counts()
    counts_df = counts.reset_index()
    counts_df.columns = ['gene', 'count']
    
    return transcripts_in_polygons, counts_df


def get_polygon_area(polygon_file):
    
    polygon_df = pd.read_csv(polygon_file)
    polygons = (
        polygon_df.groupby('polygon_id')
        .apply(lambda group: Polygon(zip(group['X'], group['Y'])))
        .reset_index(name='geometry')
    )
    polygons_gdf = gpd.GeoDataFrame(polygons, geometry='geometry')
    polygons_gdf['area'] = polygons_gdf['geometry'].area
    
    return polygons_gdf


def get_transcripts_density_in_polygon(transcripts_in_polygons, polygon_file):
    
    #transcripts_in_polygons, counts_df = get_transcripts_in_polygon(xenium_bundle, polygon_file)
    
    transcripts_in_polygons['category'] = transcripts_in_polygons['feature_name'].apply(categorize_prefix)

    counts = transcripts_in_polygons['category'].value_counts()
    
    polygons_gdf = get_polygon_area(polygon_file)
    
    return counts/polygons_gdf['area'].sum()*100


    
def get_transcripts_cor_polygon(counts_df, sudo_bulk_PATH):
    #transcripts_in_polygons, counts_df = get_transcripts_in_polygon(xenium_bundle, polygon_file)
    sudo_bulk = pd.read_csv(sudo_bulk_PATH)
    
    merged_df = counts_df.merge(sudo_bulk, how='left', on='gene')
    
    filtered_df = merged_df[
    ~merged_df['gene'].str.startswith(('UnassignedCodeword', 'NegControlCodeword', 'DeprecatedCodeword', 'NegControlProbe'))
    ]
    
    filtered_df = filtered_df[(filtered_df['count'] > 0) & (filtered_df['total_counts'] > 0)]
    corr_matrix = np.corrcoef(filtered_df['count'], filtered_df['total_counts'])
    corr = corr_matrix[0,1]
    R_sq = corr**2
    
        # Scatter plot data
    x = filtered_df['count']
    y = filtered_df['total_counts']

    # Apply log transformation
    log_x = np.log10(x)
    log_y = np.log10(y)

    # Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(log_x, log_y)

    # Regression line
    regression_line = slope * log_x + intercept

    # Plot the scatter plot and regression line
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.7)
    #plt.plot(10**log_x, 10**regression_line, color='red', label=f'Fit: y = {10**intercept:.2f}x^{slope:.2f}')
    #plt.plot(10**log_x, 10**regression_line, color='red', 
    #             label=f'Fit: y = {10**intercept:.2f}x^{slope:.2f}\nR = {corr:.3f}')
    # Annotate the R value
    plt.text(
        0.05, 0.95, f"R = {corr:.3f}", 
        transform=plt.gca().transAxes, 
        fontsize=12, 
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.5)
    )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Xenium')
    plt.ylabel('Sudo-Bulk RNA')
    plt.title('Transcript Correlation')
    #plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.show()

