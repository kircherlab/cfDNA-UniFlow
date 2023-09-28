GC uncorrected vs GC corrected cell free DNA (cfDNA) {{snakemake.wildcards.signal}} 
overlay in target region {{snakemake.wildcards.target_region}} with status: {{snakemake.wildcards.status_name}}.

## How to read the figure

The figure shows the effects of GC bias in a respective target region.

The left plot shows the uncorrected {{snakemake.wildcards.signal}} and 
the right plot the corrected {{snakemake.wildcards.signal}} around the specified target regions. 
The x-axis shows coordinates respective to the center of a target region in base pair resolution. 
The y-axis shows the average {{snakemake.wildcards.signal}} over all specified target regions. 
