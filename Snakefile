rule manuscript:
    input:
        "latex/figures/fig2-required-sample-size-comparison.pdf",
        "latex/figures/fig3-pos-prime-composition.pdf",
        "latex/figures/fig4-power-distribution-ep-approach.pdf",
        "latex/figures/fig5-power-distribution-quantile-approach.pdf",
        "latex/figures/fig6-clinical-trial-example.pdf",
        "latex/figures/fig7-matched-reward.pdf"
    output:
        "latex/main.pdf"
    shell:
        """
        set -ex
        Rscript -e "setwd('latex'); tinytex::latexmk('main.tex', bib_engine = 'biber')"
        """

rule figures:
    output:
        "latex/figures/fig2-required-sample-size-comparison.pdf",
        "latex/figures/fig3-pos-prime-composition.pdf",
        "latex/figures/fig4-power-distribution-ep-approach.pdf",
        "latex/figures/fig5-power-distribution-quantile-approach.pdf",
        "latex/figures/fig6-clinical-trial-example.pdf",
        "latex/figures/fig7-matched-reward.pdf"
    shell:
        """
        set -ex
        jupyter nbconvert \
            --to notebook \
            --execute --ExecutePreprocessor.timeout=600 \
            notebooks/figures-for-manuscript.ipynb
        """
