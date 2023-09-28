Cell free DNA (cfDNA) {{snakemake.wildcards.signal}} overlay in target region {{snakemake.wildcards.target_region}} 
comparing case: {{snakemake.wildcards.case_name}} and control: {{snakemake.wildcards.control_name}}.

## How to read the figure

The figure shows {{snakemake.wildcards.signal}} overlays of case samples in comparison to user defined controls. 
The controls can be displayed separately or as average of all controls. 
The {{snakemake.wildcards.signal}} is shown in the left plot and the corrected {{snakemake.wildcards.signal}} in the right plot.

The x-axis shows coordinates respective to the center of a target region in base pair resolution. 
The y-axis shows the average {{snakemake.wildcards.signal}} over all specified target regions.
