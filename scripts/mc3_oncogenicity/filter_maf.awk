# Filter the MC3 MAF for maf2vcf.pl. Called by run_mc3.sh as
#   awk -F'\t' -v OFS='\t' -v limit=N -v counts=FILE -v stats=FILE -f filter_maf.awk
#
# Keeps FILTER=PASS rows (the first `limit` of them when limit > 0), replaces
# every tumor/normal barcode with one placeholder pair so maf2vcf.pl writes a
# 2-sample VCF instead of one column per pair, and appends a variant_id column
# holding the GRCh37 key (chr_pos_ref_alt as the MAF writes it), which
# maf2vcf.pl copies to the VCF ID. The distinct tumors per key go to `counts`,
# row totals to `stats`.

/^#/ { next }
!hdr {
  for (i = 1; i <= NF; i++) col[$i] = i
  split("Chromosome Start_Position Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 Tumor_Sample_Barcode Matched_Norm_Sample_Barcode FILTER NCBI_Build", need, " ")
  for (k in need) if (!(need[k] in col)) { print "MAF has no " need[k] " column" > "/dev/stderr"; exit 2 }
  print $0, "variant_id"; hdr = 1; next
}
{
  rows++
  if ($col["NCBI_Build"] != "GRCh37") { print "unexpected NCBI_Build " $col["NCBI_Build"] > "/dev/stderr"; exit 2 }
  if ($col["FILTER"] != "PASS") next
  if (limit > 0 && pass >= limit) next
  pass++
  ref = $col["Reference_Allele"]; a1 = $col["Tumor_Seq_Allele1"]; a2 = $col["Tumor_Seq_Allele2"]
  alt = (a2 != ref) ? a2 : a1
  key = $col["Chromosome"] "_" $col["Start_Position"] "_" ref "_" alt
  tumor = $col["Tumor_Sample_Barcode"]
  if (!((key, tumor) in seen)) { seen[key, tumor] = 1; n[key]++ }
  if (!(key in gene)) {
    gene[key] = ("Hugo_Symbol" in col) ? $col["Hugo_Symbol"] : ""
    cls[key] = ("Variant_Classification" in col) ? $col["Variant_Classification"] : ""
    hgvsp[key] = ("HGVSp_Short" in col) ? $col["HGVSp_Short"] : ""
  }
  $col["Tumor_Sample_Barcode"] = "MC3_TUMOR"
  $col["Matched_Norm_Sample_Barcode"] = "MC3_NORMAL"
  print $0, key
}
END {
  print "key\tn_tumors\tgene\tmc3_class\thgvsp" > counts
  for (k in n) print k, n[k], gene[k], cls[k], hgvsp[k] > counts
  u = 0; for (k in n) u++
  print "maf_rows\t" rows > stats
  print "pass_rows\t" pass > stats
  print "pass_unique_variants\t" u > stats
}
