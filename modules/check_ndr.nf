process CHECK_NDR {
  input:
  path ndr_file

  output:
  val min_ndr

  script:
  """
  min_ndr=\$(awk -F',' 'NR>1 {print \$3}' ${ndr_file} | sort -n | head -n1)
  echo \$min_ndr
  """
}