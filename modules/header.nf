def logHeader(params) {
    // Log colors ANSI codes
    c_dim = "\033[2m";
    c_green = "\033[0;32m";
    c_purple = "\033[0;35m";
    c_reset = "\033[0m";

    return """-${c_dim}-------------------------------------${c_reset}-
${c_green}    ___    _   ___   _________  __ ___
   /   |  / | / / | / / ____/ |/ //   |
  / /| | /  |/ /  |/ / __/  |   // /| |
 / ___ |/ /|  / /|  / /___ /   |/ ___ |
/_/  |_/_/ |_/_/ |_/_____//_/|_/_/  |_|
                                       ${c_reset}
-${c_dim}-------------------------------------${c_reset}-
${c_purple}github.com/mlorthiois/ANNEXA${c_reset}
Reference Annotation: ${params.gtf}
Reference Genome    : ${params.fa}
Input Samplesheet   : ${params.input}
---
Filtering           : ${params.filter}
Tfkmers Model       : ${params.tfkmers_model}
Tfkmers Tokenizer   : ${params.tfkmers_tokenizer}
Tfkmers Threshold   : ${params.tfkmers_threshold}
Bambu Threshold     : ${params.bambu_threshold}
Filtering operation : ${params.operation}
-${c_dim}-------------------------------------${c_reset}-
""".stripIndent()
}
