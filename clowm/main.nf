fileEndingList = ["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz", "*.fna", "*.fna.gz", "*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]

process build_bifrost {
    tag "BuildBifrost_${name_param}"
    container 'corer:latest' // Assuming corer image contains Bifrost
    publishDir "${params.outdir}/bifrost_build", mode: "copy"

    cpus 12

    input:
      val name_param
      path inputFiles // Expects a list of files

    output:
      path "${name_param}GenomeList.txt", emit: genomelist
      path "${name_param}.gfa.gz", emit: gfa_file
      path "${name_param}.color.bfg", emit: color_file
      path "${name_param}.bfi", emit: index_file

    script:
    """
    # Create genomeList.txt with one file path per line
    printf '%s\\n' '${inputFiles.join("\n")}' > ${name_param}GenomeList.txt

    Bifrost build \
    --input-ref-file ${name_param}GenomeList.txt \
    --output-file ${name_param} \
    --colors \
    --threads ${task.cpus} \
    --kmer-length ${params.kmer_length} \
    --bloom-bits ${params.bloom_bits}  \
    ${params.clip_tips ? '--clip-tips' : ''} \
    ${params.del_isolated ? '--del-isolated' : ''} \
    """
}

process corer {
    tag "Corer_${name_param}"
    container 'corer:latest'
    publishDir "${params.outdir}/corer_analysis", mode: "copy"

    cpus 12

    input:
      val name_param
      // path original_input_files
      path gfa_file
      path color_file
      path index_file


    output:
      path "${name_param}Core.*"

    script:
    """
    Corer --igraph $gfa_file\
    --cgraph $color_file \
    --ograph ${name_param}Core \
    --threads ${task.cpus} \
    ${ params.quorum != null ? "--quorum ${params.quorum}" : ""} \
    --distance ${params.distance}
    """
}

process unzip {
  container 'python:3.13'
  input:
  path zipgenomes
  output:
  path 'output/*'

  script:
  """
  #!/usr/local/bin/python

  import zipfile
  with zipfile.ZipFile("${zipgenomes }", "r") as f:
      f.extractall("output")
  """
}

process untargz {
  container 'python:3.13'
  input:
  path zipgenomes
  output:
  path 'output/*'

  script:
  """
  #!/usr/local/bin/python

  import tarfile
  with tarfile.open("${zipgenomes}", "r") as f:
      f.extractall("output")
  """
}

workflow {
  if (params.input != null && params.bifrost_graph != null) {
    error("Specify either 'input' or 'bifrost_graph' but not both")
  }
  if (params.input != null) {
    if (params.input.endsWith(".zip")) {
    files = unzip(params.input)
    } else if (params.input.endsWith(".tar.gz")) {
    files = untargz(params.input)
    } else {
    inputChannel=Channel.fromPath(fileEndingList.collect { params.input + "/" + it },type : "file")
    files = inputChannel.collect()
    }

    bifrost_results = build_bifrost(params.dataset_name, files)

    corer(
        params.dataset_name,
        bifrost_results.gfa_file,
        bifrost_results.color_file,
        bifrost_results.index_file
    )
  } else {
    // Remove the ".gfa.gz" suffix
    String baseName = params.bifrost_graph.replace(".gfa.gz", "")

    // Generate the new filenames
    def outputBfi = "${baseName}.bfi"
    def outputColorBfg = "${baseName}.color.bfg"
    corer(
      baseName.split("/").last(),
      params.bifrost_graph,
      outputColorBfg,
      outputBfi
    )
  }

}
