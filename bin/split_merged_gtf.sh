grep "gene\." $1 | grep -P "\texon\t" > novel.gtf
grep -v "gene\." $1 | grep -P "\texon\t" > known.gtf
