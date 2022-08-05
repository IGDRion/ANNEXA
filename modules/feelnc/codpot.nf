process FEELNC_CODPOT {
  conda (params.enable_conda ? "bioconda::feelnc=0.2" : null)
  container "${ workflow.containerEngine == 'singularity' ? 
                'https://depot.galaxyproject.org/singularity/feelnc:0.2--pl526_0' : 
                'quay.io/biocontainers/feelnc:0.2--pl526_0' }"
  memory params.maxMemory

  input:
  file ref
  file fa
  file gtf

  output:
  path "feelnc_codpot_out/new.lncRNA.gtf", emit: lncRNA
  path "feelnc_codpot_out/new.mRNA.gtf", emit: mRNA

  """
  grep "protein_coding" ${ref} > known_mRNA.gtf
  grep -v "protein_coding" ${ref} > known_lncRNA.gtf

  # Source when using a container
  if [ -f "/usr/local/etc/conda/activate.d/feelnc-env.sh" ]; then
    source "/usr/local/etc/conda/activate.d/feelnc-env.sh"
  fi

  FEELnc_codpot.pl \
      -i ${gtf} \
      -g ${fa} \
      -a known_mRNA.gtf \
      -l known_lncRNA.gtf \
      --numtx=3000,3000 \
      -o new
  """
}

