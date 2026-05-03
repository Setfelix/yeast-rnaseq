process INFER_STRANDEDNESS {
  tag "$sample"

  publishDir "${params.outdir}/qc/strandedness", mode: 'copy'

  input:
  tuple val(sample), path(bam), path(bai), val(strandedness), val(condition)
  path annotation_bed12

  output:
  tuple val(sample), path(bam), path(bai), val(condition), path("${sample}.strandedness.txt"), emit: resolved
  path "${sample}.infer_experiment.txt", emit: report

  script:
  """
  case "${strandedness}" in
    forward|reverse|unstranded)
      printf "%s\n" "${strandedness}" > "${sample}.strandedness.txt"
      printf "strandedness_source: provided\nstrandedness: %s\n" "${strandedness}" > "${sample}.infer_experiment.txt"
      ;;
    ""|auto)
      infer_experiment.py \
        -i "${bam}" \
        -r "${annotation_bed12}" \
        > "${sample}.infer_experiment.txt"

      python3 - "${sample}.infer_experiment.txt" "${params.strandedness_inference_threshold}" "${params.strandedness_unstranded_tolerance}" > "${sample}.strandedness.txt" <<'PY'
import pathlib
import re
import sys

report_path = pathlib.Path(sys.argv[1])
threshold = float(sys.argv[2])
tolerance = float(sys.argv[3])
text = report_path.read_text()

forward_match = re.search(r'Fraction of reads explained by "1\\+\\+,1--,2\\+-,2-\\+":\\s*([0-9.]+)', text)
reverse_match = re.search(r'Fraction of reads explained by "1\\+-,1-\\+,2\\+\\+,2--":\\s*([0-9.]+)', text)

if not forward_match or not reverse_match:
    raise SystemExit(f"Could not parse infer_experiment.py output from {report_path}")

forward = float(forward_match.group(1))
reverse = float(reverse_match.group(1))

if forward >= threshold and forward > reverse:
    print("forward")
elif reverse >= threshold and reverse > forward:
    print("reverse")
elif abs(forward - reverse) <= tolerance:
    print("unstranded")
else:
    raise SystemExit(
        f"Ambiguous strandedness inference for {report_path.name}: "
        f"forward={forward:.4f}, reverse={reverse:.4f}. "
        f"Provide strandedness explicitly or adjust the inference thresholds."
    )
PY
      ;;
    *)
      echo "Unsupported strandedness '${strandedness}'. Use one of: unstranded, forward, reverse" >&2
      exit 1
      ;;
  esac
  """

  stub:
  """
  if [ -n "${strandedness}" ]; then
    resolved="${strandedness}"
    printf "strandedness_source: provided\nstrandedness: %s\n" "${strandedness}" > "${sample}.infer_experiment.txt"
  else
    resolved="reverse"
    printf 'This is PairEnd Data\nFraction of reads failed to determine: 0.02\nFraction of reads explained by "1++,1--,2+-,2-+": 0.03\nFraction of reads explained by "1+-,1-+,2++,2--": 0.95\n' > "${sample}.infer_experiment.txt"
  fi

  printf "%s\n" "\$resolved" > "${sample}.strandedness.txt"
  """
}
