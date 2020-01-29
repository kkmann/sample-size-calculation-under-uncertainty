rule introduction:
    input:
        notebook = "notebooks/required-sample-size-vs-prior.ipynb"
    output:
        "latex/figures/power-constraint-comparison.pdf"
    shell:
        """
        jupyter nbconvert --execute {input.notebook}
        rm notebooks/*.html
        """
