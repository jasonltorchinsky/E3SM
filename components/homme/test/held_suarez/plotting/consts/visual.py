# Standard Library Imports

# Third-Party Library Imports

# Local Library Imports

diff_cmap: str = "RdBu_r" # Colormap for difference quantities
heating_cmap: str = "hot" # Colormap for heating-related quanities
flux_cmap: str = "magma" # Colormap for flux-related quantities
lw_cmap: str = "Blues" # Colormap for liquid water-related values
iw_cmap: str = "Purples" # Colormap for ice-water related quantities
cw_cmap: str = "Greens" # Colormap for cloud water-related quantities
plot_colors: list[str] = ["#332288", "#117733", "#44AA99", "#88CCEE", 
    "#DDCC77", "#CC6677", "#AA4499", "#882255"] # Color-blind-friendly color palette
plot_markers: list[str] = ["o", "^", "s", "D", "p", "H", "P", "*"] # Marker styles for scatter plots