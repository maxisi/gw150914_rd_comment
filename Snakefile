rule finchmoore
    output:
        "src/data/fm_t0_kde_grid.h5"
    script:
        "src/scripts/finchmoore.py"
