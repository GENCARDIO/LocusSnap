#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/locus_snap_matplotlib"
mkdir -p "$MPLCONFIGDIR"

python3 generate_demo_data.py

python3 -m locus_snap \
  --bam test/test.bam \
  --sample_label 'Primary locus with inferred discordant-mate window' \
  --region chr9:101867481-101867620 \
  --mate_view \
  --mate_window_source discordant \
  --layout expand \
  --display_mode squish \
  --max_alignment_depth 0 \
  --refseq none \
  --output_dir out \
  --output_name 06_mate_view_discordant \
  --fig_width 14 \
  --dpi 120

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --custom_track 'out/demo_data/annotations/demo_regions.bed,bed,Candidate regions,#000000,pack,0.42' \
  --custom_track 'out/demo_data/annotations/demo_genes.gtf,gtf,GENCODE genes,#17217a,pack,0.72' \
  --custom_track 'out/demo_data/variants/demo_variants.vcf.gz,vcf,Expanded variants,#7a1f5c,pack,0.42' \
  --display_mode collapse \
  --refseq none \
  --output_dir out \
  --output_name 11_custom_track_definitions \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam test/test.bam \
  --sample_label 'Visible mates linked on shared rows' \
  --region chr9:101867481-101867620 \
  --view_as_pairs \
  --layout pack \
  --display_mode expand \
  --max_alignment_depth 0 \
  --refseq none \
  --output_dir out \
  --output_name 15_view_as_pairs \
  --fig_width 14 \
  --dpi 120

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --display_mode squish \
  --layout pack \
  --coverage_vaf_threshold 0.20 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 16_coverage_snv_vaf \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_cnv_tumour.bam \
  --sample_label 'Tumour · 75% purity · 16,000 reads · CN1 loss and CN3 gain' \
  --fasta out/demo_data/reference/demo_cnv_reference.fa \
  --region chrCNV:1-80000 \
  --custom_track 'out/demo_data/annotations/demo_cnv_states.bed,bed,Copy-number state,#333333,pack,0.42' \
  --track out/demo_data/annotations/demo_cnv.seg \
  --track_label 'Tumour CNV · purity-adjusted log2 ratio' \
  --baf_vcf out/demo_data/variants/demo_baf.vcf.gz \
  --baf_sample Tumour \
  --baf_track_label 'Tumour BAF · CN1 loss: 0.20/0.80 · CN3 gain: 0.36/0.64' \
  --display_mode squish \
  --layout pack \
  --no_alignments \
  --grid_mode bands \
  --max_alignment_depth 0 \
  --coverage_vaf_threshold 0.12 \
  --min_baseq 20 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 18_variant_evidence_baf_loh \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --haplotype_view split \
  --display_mode squish \
  --layout pack \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 19_haplotype_split_view \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --sample_label 'Tumour reads grouped by nucleotide at the selected SNV' \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --layout expand \
  --display_mode expand \
  --sort_by base \
  --sort_base_position 119 \
  --max_rows 48 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 22_sort_by_snv_base \
  --fig_width 13 \
  --dpi 140

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_met_ex14.bam \
  --sample_label 'METex14-positive lung adenocarcinoma · synthetic RNA-seq' \
  --region chr7:116771401-116775200 \
  --custom_track 'out/demo_data/annotations/demo_met_ex14.gtf,gtf,MET · NM_000245.4 (exons 13–15),#17217a,collapse,0.52' \
  --custom_track 'out/demo_data/variants/demo_met_ex14.vcf.gz,vcf,MET c.3028+1G>T · exon 14 donor,#7a1f5c,pack,0.42' \
  --sashimi \
  --min_junction_reads 5 \
  --sashimi_strand combined \
  --display_mode squish \
  --layout pack \
  --max_alignment_depth 120 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 24_rnaseq_sashimi \
  --fig_width 14 \
  --dpi 140

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101865501-101869500 \
  --config out/demo_data/config/demo_chipseq.yaml \
  --sample_label 'CTCF ChIP-seq normalized signal' \
  --custom_track 'out/demo_data/signals/demo_ctcf_control.signal.gz,signal,Wehi-CT control,#00695c,collapse,0.95' \
  --custom_track 'out/demo_data/signals/demo_ctcf_knockdown.signal.gz,signal,Wehi-TFII-I-KD,#22d3a6,collapse,0.95' \
  --custom_track 'out/demo_data/signals/demo_ctcf_mel.signal.gz,signal,MEL CTCF,#4d9bd6,collapse,0.95' \
  --custom_track 'out/demo_data/annotations/demo_ctcf_genes.gtf,gtf,CTCF target genes,#17217a,collapse,0.68' \
  --display_mode collapse \
  --no_alignments \
  --no_coverage \
  --refseq none \
  --output_dir out \
  --output_name 26_chipseq_peaks_density \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --bam out/demo_data/alignments/demo_normal.bam \
  --bam out/demo_data/alignments/demo_relapse.bam \
  --sample_label Tumour \
  --sample_label Normal \
  --sample_label Relapse \
  --vcf_companion out/demo_data/variants/demo_tumour.vcf.gz \
  --vcf_companion none \
  --vcf_companion out/demo_data/variants/demo_relapse.vcf.gz \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:81-180 \
  --display_mode squish \
  --layout pack \
  --max_alignment_depth 40 \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 27_multi_bam_vcf_companions \
  --fig_width 14 \
  --dpi 150

python3 -m locus_snap \
  --bam test/test.bam \
  --region chr9:101867481-101867620 \
  --output_dir out \
  --output_name 30_default_refseq_isoforms \
  --fig_width 14 \
  --dpi 100

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_structural_variants.bam \
  --sample_label 'Tumour · reciprocal translocation evidence' \
  --region chr1:7351-7850 \
  --region chr2:4751-5250 \
  --region_label 'chr1 breakpoint · BND partner chr2' \
  --region_label 'chr2 breakpoint · BND partner chr1' \
  --link_breakpoints \
  --custom_track 'out/demo_data/variants/demo_structural_variants.vcf.gz,vcf,Somatic BND,#7a1f5c,collapse,0.34' \
  --only discordant split softclip \
  --min_softclip 20 \
  --view_as_pairs \
  --layout pack \
  --display_mode squish \
  --sort_by gap_length \
  --sort_order desc \
  --max_rows 40 \
  --max_alignment_depth 100 \
  --no_annotate \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 34_explicit_multilocus_breakpoint \
  --fig_width 16 \
  --dpi 150

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_structural_variants.bam \
  --sample_label 'Tumour · deletion, tandem duplication, inversion, and chr1–chr2 translocation' \
  --region chr1:1001-8500 \
  --custom_track 'out/demo_data/variants/demo_structural_variants.vcf.gz,vcf,Somatic SVs · DEL · DUP · INV · TRA,#7a1f5c,pack,0.46' \
  --min_softclip 20 \
  --view_as_pairs \
  --layout pack \
  --display_mode squish \
  --sort_by gap_length \
  --sort_order desc \
  --max_rows 80 \
  --max_alignment_depth 120 \
  --no_annotate \
  --genome none \
  --refseq none \
  --output_dir out \
  --output_name 31_structural_variant_evidence \
  --fig_width 16 \
  --dpi 150

python3 -m locus_snap \
  --bam test/test.bam \
  --sample_label 'Close zoom · soft-clipped sequence shown base by base' \
  --region chr9:101867541-101867580 \
  --config out/demo_data/config/demo_softclip_zoom.yaml \
  --only softclip \
  --min_softclip 1 \
  --layout expand \
  --display_mode expand \
  --sort_by start \
  --sort_order asc \
  --max_alignment_depth 0 \
  --no_coverage \
  --no_ideogram \
  --refseq none \
  --output_dir out \
  --output_name 32_close_zoom_softclipped_bases \
  --fig_width 10 \
  --dpi 150

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_insertions.bam \
  --sample_label 'Close zoom · short insertions marked at their breakpoint' \
  --fasta out/demo_data/reference/demo_reference.fa \
  --region chrDemo:101-140 \
  --config out/demo_data/config/demo_softclip_zoom.yaml \
  --only gapped \
  --layout expand \
  --display_mode expand \
  --sort_by start \
  --sort_order asc \
  --max_alignment_depth 0 \
  --no_annotate \
  --no_coverage \
  --no_ideogram \
  --refseq none \
  --output_dir out \
  --output_name 33_close_zoom_insertions \
  --fig_width 10 \
  --dpi 150

python3 -m locus_snap \
  --bam out/demo_data/alignments/demo_tumour.bam \
  --bam out/demo_data/alignments/demo_normal.bam \
  --bam out/demo_data/alignments/demo_relapse.bam \
  --sample_label Tumour \
  --sample_label Normal \
  --sample_label Relapse \
  --vcf_companion out/demo_data/variants/demo_tumour.vcf.gz \
  --vcf_companion none \
  --vcf_companion out/demo_data/variants/demo_relapse.vcf.gz \
  --fasta out/demo_data/reference/demo_reference.fa \
  --batch_regions out/demo_data/annotations/demo_multi_sample_review.bed \
  --report 35_multi_sample_batch_report.html \
  --threads 3 \
  --display_mode squish \
  --layout pack \
  --max_alignment_depth 40 \
  --genome none \
  --refseq none \
  --output_dir out/multi_sample_batch \
  --fig_width 14 \
  --dpi 120

if command -v google-chrome >/dev/null 2>&1; then
  google-chrome \
    --headless \
    --disable-gpu \
    --no-sandbox \
    --hide-scrollbars \
    --window-size=1600,1100 \
    --screenshot="$PROJECT_DIR/out/35_multi_sample_batch_report.png" \
    "file://$PROJECT_DIR/out/multi_sample_batch/35_multi_sample_batch_report.html"
else
  printf '%s\n' 'Skipped report preview: google-chrome is not installed.'
fi

printf '%s\n' 'Regenerated the curated README figures and multi-sample batch report in out/'
