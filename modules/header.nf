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
${c_purple}github.com/igdrion/ANNEXA${c_reset}
ANNEXA version        : ${workflow.manifest.version}
---
Reference Annotation  : ${params.gtf}
Reference Genome      : ${params.fa}
Input Samplesheet     : ${params.input}
---
Transcript discovery  : ${params.tx_discovery}
Filtering             : ${params.filter}
Tfkmers Model         : ${params.tfkmers_model}
Tfkmers Tokenizer     : ${params.tfkmers_tokenizer}
Tfkmers Threshold     : ${params.tfkmers_threshold}
Bambu Threshold       : ${params.bambu_threshold}
Bambu Single Exons    : ${params.bambu_singleexon}
Bambu Recommended NDR : ${params.bambu_rec_ndr}
Filtering operation   : ${params.operation}
Stranded              : ${params.bambu_strand}
-${c_dim}-------------------------------------${c_reset}-
""".stripIndent()
}

def helpHeader() {
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
${c_purple}github.com/igdrion/ANNEXA${c_reset}
ANNEXA version : ${workflow.manifest.version}
""".stripIndent()
}