def plot_samples_with_grid(adata, samples, group, concatenated_sdata):
    import matplotlib.pyplot as plt
    import numpy as np

    for sample in samples:
        shapes = concatenated_sdata[f'{sample}_cell_boundaries'].copy()
        mask = adata.obs['sample'] == sample
        obs_subset = adata.obs[mask][['instance_id'] + group]
        shapes = shapes.merge(obs_subset, left_index=True, right_on='instance_id', how='inner')

        x_min, x_max = shapes.geometry.centroid.x.min(), shapes.geometry.centroid.x.max()
        y_min, y_max = shapes.geometry.centroid.y.min(), shapes.geometry.centroid.y.max()

        for i in group:
            fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
            shapes.plot(column=i, ax=ax, legend=True, cmap='tab20', linewidth=0,
                        legend_kwds={'loc': 'upper left', 'bbox_to_anchor': (1, 1)})

            # major grid every 500 units
            major_grid_spacing = 500
            x_major_ticks = np.arange(0, x_max + major_grid_spacing, major_grid_spacing)
            y_major_ticks = np.arange(0, y_max + major_grid_spacing, major_grid_spacing)

            # minor grid every 100 units
            minor_grid_spacing = 100
            x_minor_ticks = np.arange(0, x_max + minor_grid_spacing, minor_grid_spacing)
            y_minor_ticks = np.arange(0, y_max + minor_grid_spacing, minor_grid_spacing)

            ax.set_xticks(x_major_ticks)
            ax.set_yticks(y_major_ticks)
            ax.set_xticks(x_minor_ticks, minor=True)
            ax.set_yticks(y_minor_ticks, minor=True)

            ax.grid(which='major', color='red', linestyle='-', linewidth=1.5, alpha=0.7)
            ax.grid(which='minor', color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

            for x_pos in x_major_ticks[::2]:
                ax.text(x_pos, y_max + 50, f'{int(x_pos)}',
                        ha='center', va='bottom', fontsize=10, color='red', fontweight='bold')
            for y_pos in y_major_ticks[::2]:
                ax.text(x_min - 50, y_pos, f'{int(y_pos)}',
                        ha='right', va='center', fontsize=10, color='red', fontweight='bold')

            ax.set_title(f'{sample} - {i}')
            ax.invert_yaxis()
            plt.tight_layout()
            plt.show()


def crop_with_grids(adata):
    import matplotlib.pyplot as plt
    import numpy as np

    labels = adata.obs['leiden_harmony_0.7'].astype(str)
    x = adata.obs['array_col'].values
    y = adata.obs['array_row'].values

    unique_labels = sorted(labels.unique(), key=lambda x: int(x))
    colors = plt.cm.tab20(range(len(unique_labels)))
    label_to_color = dict(zip(unique_labels, colors))

    # create figure with grid
    fig, ax = plt.subplots(figsize=(14, 14))

    # plot the data
    for label in unique_labels:
        mask = labels == label
        ax.scatter(x[mask], y[mask], color=label_to_color[label], s=0.5, label=label, alpha=0.8)

    # add grid lines every 250 units
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    # major grid every 500 units
    major_grid_spacing = 500
    x_major_ticks = np.arange(0, x_max + major_grid_spacing, major_grid_spacing)
    y_major_ticks = np.arange(0, y_max + major_grid_spacing, major_grid_spacing)

    # minor grid every 100 units
    minor_grid_spacing = 100
    x_minor_ticks = np.arange(0, x_max + minor_grid_spacing, minor_grid_spacing)
    y_minor_ticks = np.arange(0, y_max + minor_grid_spacing, minor_grid_spacing)

    # set ticks
    ax.set_xticks(x_major_ticks)
    ax.set_yticks(y_major_ticks)
    ax.set_xticks(x_minor_ticks, minor=True)
    ax.set_yticks(y_minor_ticks, minor=True)

    # add grid
    ax.grid(which='major', color='red', linestyle='-', linewidth=1.5, alpha=0.7)
    ax.grid(which='minor', color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

    # add coordinate labels on grid lines
    for x_pos in x_major_ticks[::2]:  # label every other major tick to avoid crowding
        ax.text(x_pos, y_max + 50, f'{int(x_pos)}',
                ha='center', va='bottom', fontsize=10, color='red', fontweight='bold')

    for y_pos in y_major_ticks[::2]:
        ax.text(x_min - 50, y_pos, f'{int(y_pos)}',
                ha='right', va='center', fontsize=10, color='red', fontweight='bold')

    # add corner coordinates for reference
    ax.text(x_min, y_min, f'({int(x_min)}, {int(y_min)})',
            ha='left', va='bottom', fontsize=12, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(x_max, y_min, f'({int(x_max)}, {int(y_min)})',
            ha='right', va='bottom', fontsize=12, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(x_min, y_max, f'({int(x_min)}, {int(y_max)})',
            ha='left', va='top', fontsize=12, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.text(x_max, y_max, f'({int(x_max)}, {int(y_max)})',
            ha='right', va='top', fontsize=12, color='blue', fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # # add legend
    # ax.legend(markerscale=5, bbox_to_anchor=(1.05, 1),
    #             loc='upper left', title='Cluster', fontsize=8)

    ax.set_aspect('equal')
    ax.set_xlabel('array_col', fontsize=12, fontweight='bold')
    ax.set_ylabel('array_row', fontsize=12, fontweight='bold')
    # ax.set_title('BANKSY Spatial Clusters with Coordinate Grid\n(Red lines: every 500 units, Gray dots: every 100 units)',
    #              fontsize=14, fontweight='bold')

    plt.tight_layout()
    # plt.savefig('test_with_grid.png', dpi=300, bbox_inches='tight')
    # print(f"Grid visualization saved as 'test_with_grid.png'")
    print(f"\nCurrent data range:")
    print(f"  array_col (x): {x_min:.0f} to {x_max:.0f}")
    print(f"  array_row (y): {y_min:.0f} to {y_max:.0f}")
    print(f"\nRed grid lines mark every {major_grid_spacing} units")
    print(f"Gray dotted lines mark every {minor_grid_spacing} units")
    plt.show()
