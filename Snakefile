rule finchmoore1:
    output:
        "src/data/A1_posterior_tref.dat"
    script:
        "src/scripts/finchmoore1.py"

rule finchmoore2:
    input:
        "src/data/A1_posterior_tref.dat"
    output:
        "src/data/fm_t0_kde_grid.h5"
    script:
        "src/scripts/finchmoore2.py"
