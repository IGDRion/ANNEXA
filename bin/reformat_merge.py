#!/usr/bin/env python3

tracking_file = "gffcmp.tracking"
gtf_file = "gffcmp.combined.gtf"
output_file = "reformated.gtf"

t_map = {}     # transcript replacement
g_map = {}     # gene replacement
extra_map = {} # extra attributes (validated, alt_transcript)

# --- Step 1: Parse tracking file ---
with open(tracking_file) as tf:
    for line in tf:
        if not line.strip():
            continue
        fields = line.strip().split("\t")
        tcons = fields[0]   # e.g. TCONS_xxx
        xloc  = fields[1]   # e.g. XLOC_xxx

        q1, q2 = None, None
        for f in fields:
            if f.startswith("q1:"):
                q1 = f
            elif f.startswith("q2:"):
                q2 = f

        gene_id, tx_id = None, None
        extras = []

        if q1 and q1 != "q1:-":
            parts1 = q1.split(":", 1)[1].split("|")
            if len(parts1) >= 2:
                gene_id, tx_id = parts1[0], parts1[1]

            validate_tag = None

            if q2 and q2 != "q2:-":
                parts2 = q2.split(":", 1)[1].split("|")
                if len(parts2) >= 2:
                    gene2, tx2 = parts2[0], parts2[1]
                    if not (gene2 == gene_id and tx2 == tx_id):
                        extras.append('validated "both"')
                        extras.append(f'alt_transcript "{tx2}"')
                else:
                    tx2 = None
                # No validate tag if q2 present, only validated "both" if differ
            else:
                # Only q1: decide validate tag based on transcript prefix
                if tx_id.startswith("BambuTx"):
                    validate_tag = 'validate "bambu"'
                elif tx_id.startswith("MSTRG."):
                    validate_tag = 'validate "stringtie"'

            if validate_tag:
                extras.append(validate_tag)

        elif q2 and q2 != "q2:-":
            parts2 = q2.split(":", 1)[1].split("|")
            if len(parts2) >= 2:
                gene_id, tx_id = parts2[0], parts2[1]
                # q2 alone: no validate tag (per your request)

        if gene_id and tx_id:
            t_map[tcons] = tx_id
            g_map[xloc]  = gene_id
            if extras:
                extra_map[(tcons, xloc)] = "; ".join(extras)

# --- Step 2: Rewrite GTF ---
with open(gtf_file) as gf, open(output_file, "w") as out:
    for line in gf:
        if line.startswith("#") or not line.strip():
            out.write(line)
            continue

        orig_line = line.rstrip("\n")
        fields = orig_line.split("\t")
        if len(fields) < 9:
            out.write(line)
            continue

        attrs = fields[8]

        # Replace transcript_id
        for old_t, new_t in t_map.items():
            if f'transcript_id "{old_t}"' in attrs:
                attrs = attrs.replace(f'transcript_id "{old_t}"',
                                      f'transcript_id "{new_t}"')

        # Replace gene_id
        for old_g, new_g in g_map.items():
            if f'gene_id "{old_g}"' in attrs:
                attrs = attrs.replace(f'gene_id "{old_g}"',
                                      f'gene_id "{new_g}"')

        # Remove unwanted attributes
        attr_parts = [a.strip() for a in attrs.split(";") if a.strip()]
        new_attrs = []
        for a in attr_parts:
            if (a.startswith("oId ") or 
                a.startswith("tss_id ") or 
                a.startswith("num_samples ") or 
                a.startswith("contained_in ") or
                a == 'validated "q2"'):  # remove validated "q2"
                continue
            new_attrs.append(a)

        attrs = "; ".join(new_attrs)

        # Add extras if needed (key off original TCONS/XLOC)
        t_id = None
        g_id = None
        if 'transcript_id "' in line:
            t_id = line.split('transcript_id "')[1].split('"')[0]
        if 'gene_id "' in line:
            g_id = line.split('gene_id "')[1].split('"')[0]
        if t_id and g_id and (t_id, g_id) in extra_map:
            attrs = attrs + "; " + extra_map[(t_id, g_id)]

        # Reorder: transcript_id first, then gene_id, then others
        attr_list = [a.strip() for a in attrs.split(";") if a.strip()]
        ordered = []
        for key in ["transcript_id", "gene_id"]:
            selected = [x for x in attr_list if x.startswith(key)]
            ordered.extend(selected)
            attr_list = [x for x in attr_list if not x.startswith(key)]
        ordered.extend(attr_list)

        fields[8] = "; ".join(ordered) + ";"
        out.write("\t".join(fields) + "\n")
