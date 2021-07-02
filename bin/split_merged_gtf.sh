grep "gene\." $1 | awk '$3=="exon"' > novel.gtf
grep -v "gene\." $1 | awk '$3=="exon"' > known.gtf
