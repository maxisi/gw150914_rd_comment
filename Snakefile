rule fig1prep:
    input:
        "src/data/figure1/figure1/pyring_t0-bad_timestamps-bad_seglen-bad/result.nc",
        "src/data/figure1/figure1/pyring_t0-good_timestamps-good_seglen-good/result.nc",
        "src/data/figure1/figure1/ringdown_t0-bad_timestamps-good_seglen-bad/result.nc",
        "src/data/figure1/figure1/ringdown_t0-bad_timestamps-good_seglen-good/result.nc",
        "src/data/figure1/figure1/ringdown_t0-good_timestamps-good_seglen-good/result.nc"
    output:
        "src/data/fig1_samples_dict.pkl"
    cache:
        True
    script:
        "src/scripts/figure1_prep.py"

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
