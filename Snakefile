rule finchmoore1:
    input:
        "src/data/finch_moore/dataset/posterior_samples/GW150914/3W220221/posterior_samples.dat"
    output:
        "src/data/A1_posterior_tref.dat"
    cache:
        True
    script:
        "src/scripts/finchmoore1.py"

rule finchmoore2:
    input:
        "src/data/A1_posterior_tref.dat"
    output:
        "src/data/fm_t0_kde_grid.h5"
    cache:
        True
    script:
        "src/scripts/finchmoore2.py"
