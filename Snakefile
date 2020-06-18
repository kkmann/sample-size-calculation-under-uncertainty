rule figures:
    output:
        "latex/figures/power-constraint-comparison.pdf",
        "latex/figures/pos-components.pdf",
        "latex/figures/quantile-vs-ep.pdf",
        "latex/figures/power-distribution-examples.pdf",
        "latex/figures/power-distribution-quantile-approach.pdf",
        "latex/figures/matched-reward.pdf"
    shell:
        """
        set -ex
        mkdir -p latex/figures

        Rscript R/figure-2_required-sample-size-vs-prior.R
        Rscript R/figure-3_pos-composition.R
        Rscript R/figure-4_distribution-of-power-ep-approach.R
        Rscript R/figure-5-and-6_distribution-of-power-quantile-approach.R
        Rscript R/figure-7_utility-matching.R
        """
