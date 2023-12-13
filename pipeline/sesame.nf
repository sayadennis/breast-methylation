params.idat_dir = "/projects/p30791/methylation/data/IDAT_all"
params.meta_fn = "/projects/p30791/methylation/data/meta.csv"
params.out_dir = "/projects/p30791/methylation/sesame_out"
params.plot_dir = "/projects/p30791/methylation/plots"

process summarize_meta {
    cpus 1
    memory '3 GB'
    time '1h'
    queue 'short'

    output:
    stdout

    script:
    """
    module purge all
    module load python-miniconda3/4.12.0
    source activate methylation

    python $HOME/breast-methylation/pipeline/00_summarize_meta.py ${params.meta_fn} ${params.out_dir}/data_summary
    """
}

process preprocess {

    input:
    val result

    output:
    stdout

    script:
    """
    module purge all
    module load R/4.3.0

    echo "Welcome to the SeSaMe Nextflow pipeline!"
    echo "The IDAT directory is ${params.idat_dir}"
    echo "The output directory is ${params.out_dir}"
    Rscript --vanilla $HOME/breast-methylation/pipeline/01_preprocess_qc.R ${params.idat_dir} ${params.out_dir}
    """
}

process plot_qc {
    cpus 1
    memory '48 GB'
    time '24h'
    queue 'normal'

    input:
    val result

    output:
    stdout

    script:
    """
    module purge all
    module load R/4.3.0

    echo "Plotting QC metrics... Resulting plots will be saved to ${params.plot_dir}"
    Rscript --vanilla $HOME/breast-methylation/pipeline/02_plot_quality.R ${params.out_dir} ${params.plot_dir} ${params.meta_fn}
    """
}

process infer_meta {
    cpus 1
    memory '96 GB'
    time '4h'
    queue 'short'

    input:
    val result

    script:
    """
    module purge all
    module load R/4.3.0

    Rscript --vanilla $HOME/breast-methylation/pipeline/03_meta_inference.R ${params.out_dir} ${params.meta_fn} ${params.out_dir}/data_summary
    """
}

workflow {
    summarize_meta | preprocess | plot_qc | infer_meta
}

