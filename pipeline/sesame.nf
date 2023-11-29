params.idat_dir = "/projects/p30791/methylation/copied_from_b1122/data/IDAT_Files/IDAT_only"
params.out_dir = "/projects/p30791/methylation/sesame_out"

process basic_preprocess {

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

workflow {
  basic_preprocess | view
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary 
        --------------------------- 
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: "sayarenedennis@northwestern.edu", subject: "nextflow-sesame-job", body: msg)
}
