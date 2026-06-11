from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from massflow.runtime_config import load_runtime_config

from ._utils import make_json_safe


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="massflow.agent_tools",
        description="Local debug entry for the first MassFlow agent tools.",
    )
    parser.add_argument(
        "--config",
        help="Optional path to massflow.toml. When omitted, the current directory is searched.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    load_parser = subparsers.add_parser("load", help="Load imzML header and dataset index.")
    load_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    load_parser.add_argument("--output-dir", help="Directory used to store the JSON summary artifact.")
    load_parser.add_argument("--max-threads", type=int, help="Maximum worker threads for MassFlow.")
    _add_biological_metadata_args(load_parser)

    qc_parser = subparsers.add_parser("qc", help="Run lightweight QC on sampled spectra.")
    qc_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    qc_parser.add_argument("--output-dir", help="Directory used to store the JSON summary artifact.")
    qc_parser.add_argument("--max-threads", type=int, help="Maximum worker threads for MassFlow.")
    qc_parser.add_argument("--sample-size", type=int, help="Number of sampled spectra for QC.")
    qc_parser.add_argument(
        "--noise-method",
        help="Noise estimation method forwarded to SpectrumPreprocess.",
    )
    _add_biological_metadata_args(qc_parser)

    recommend_parser = subparsers.add_parser(
        "recommend-preprocess",
        help="Recommend preprocessing operations from dataset QC evidence.",
    )
    recommend_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    recommend_parser.add_argument("--output-dir", help="Directory used to store the JSON summary artifact.")
    recommend_parser.add_argument("--max-threads", type=int, help="Maximum worker threads for MassFlow.")
    recommend_parser.add_argument("--sample-size", type=int, help="Number of sampled spectra for QC.")
    recommend_parser.add_argument(
        "--noise-method",
        help="Noise estimation method forwarded to SpectrumPreprocess.",
    )
    _add_biological_metadata_args(recommend_parser)

    preprocess_parser = subparsers.add_parser("preprocess", help="Run preprocessing pipeline.")
    preprocess_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    preprocess_parser.add_argument("--output-dir", help="Directory used to store processed outputs.")
    preprocess_parser.add_argument("--batch-size", type=int, help="Batch size for flat pipeline.")
    preprocess_parser.add_argument("--max-threads", type=int, help="Maximum worker threads for MassFlow.")
    preprocess_parser.add_argument("--qc-sample-size", type=int, help="QC sample size before/after preprocessing.")
    preprocess_parser.add_argument(
        "--operations-json",
        help="JSON list of preprocessing operations. Example: "
        "'[{\"name\":\"noise_reduction\",\"params\":{\"method\":\"gaussian_numba\"}}]'",
    )

    spatial_parser = subparsers.add_parser("spatial-preview", help="Generate TIC and base-peak spatial preview images.")
    spatial_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    spatial_parser.add_argument("--output-dir", help="Directory used to store spatial preview outputs.")
    spatial_parser.add_argument("--max-threads", type=int, help="Maximum worker threads for MassFlow.")
    spatial_parser.add_argument(
        "--colormap",
        default="magma",
        help="Matplotlib colormap used when saving preview images.",
    )

    mask_parser = subparsers.add_parser("apply-mask", help="Apply a PNG mask to imzML coordinates.")
    mask_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    mask_parser.add_argument("--mask-png", required=True, help="Path to the PNG mask file.")
    mask_parser.add_argument("--output-dir", help="Directory used to store mask application outputs.")
    mask_parser.add_argument("--max-threads", type=int, help="Maximum worker threads for MassFlow.")
    mask_parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Threshold used to binarize the mask image.",
    )
    mask_parser.add_argument(
        "--keep-mode",
        choices=["white_foreground", "black_foreground"],
        default="white_foreground",
        help="Interpretation of the foreground pixels in the PNG mask.",
    )

    workflow_parser = subparsers.add_parser("workflow", help="Run the end-to-end bootstrap workflow.")
    workflow_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    workflow_parser.add_argument("--output-dir", help="Directory used to store workflow outputs.")
    workflow_parser.add_argument("--no-preprocess", action="store_true", help="Skip preprocessing stage.")
    workflow_parser.add_argument(
        "--operations-json",
        help="JSON list of preprocessing operations for the workflow stage.",
    )
    _add_biological_metadata_args(workflow_parser)
    _add_chemical_annotation_args(workflow_parser)

    validate_parser = subparsers.add_parser("validate", help="Compare integrated conclusion with a legacy statistical baseline.")
    validate_parser.add_argument("--input", required=True, help="Path to the target .imzML file.")
    validate_parser.add_argument("--output-dir", help="Directory used to store validation outputs.")
    validate_parser.add_argument("--no-preprocess", action="store_true", help="Skip preprocessing stage.")
    validate_parser.add_argument(
        "--operations-json",
        help="JSON list of preprocessing operations for the validation workflow stage.",
    )
    _add_biological_metadata_args(validate_parser)
    _add_chemical_annotation_args(validate_parser)

    doctor_parser = subparsers.add_parser("doctor", help="Run environment and input self-checks.")
    doctor_parser.add_argument("--input", help="Optional path to the target .imzML file.")
    doctor_parser.add_argument("--output-dir", help="Optional writable output directory to probe.")
    doctor_parser.add_argument("--require-mcp", action="store_true", help="Fail when optional MCP dependency is missing.")

    usage_parser = subparsers.add_parser("usage-modes", help="Print the recommended switching method for the two MCP usage modes.")
    usage_parser.add_argument("--json", action="store_true", help="Print usage instructions as JSON.")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    config = load_runtime_config(args.config)

    try:
        if args.command == "load":
            from .qc_tools import load_imzml

            result = load_imzml(
                args.input,
                biological_metadata_json=args.biological_metadata_json,
                biological_metadata_path=args.biological_metadata_path,
                max_threads=args.max_threads or config.agent_tools.max_threads,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        elif args.command == "qc":
            from .qc_tools import dataset_qc_summary

            result = dataset_qc_summary(
                args.input,
                biological_metadata_json=args.biological_metadata_json,
                biological_metadata_path=args.biological_metadata_path,
                max_threads=args.max_threads or config.agent_tools.max_threads,
                sample_size=args.sample_size or config.agent_tools.qc_sample_size,
                noise_method=args.noise_method or config.agent_tools.noise_method,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        elif args.command == "preprocess":
            from .preprocess_tools import run_preprocessing

            operations = json.loads(args.operations_json) if args.operations_json else None
            result = run_preprocessing(
                args.input,
                batch_size=args.batch_size or config.agent_tools.batch_size,
                max_threads=args.max_threads or config.agent_tools.max_threads,
                qc_sample_size=args.qc_sample_size or config.agent_tools.preprocess_qc_sample_size,
                operations=operations,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        elif args.command == "spatial-preview":
            from .spatial_tools import build_spatial_preview

            result = build_spatial_preview(
                args.input,
                max_threads=args.max_threads or config.agent_tools.max_threads,
                colormap=args.colormap,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        elif args.command == "apply-mask":
            from .mask_tools import apply_spatial_png_mask

            result = apply_spatial_png_mask(
                args.input,
                args.mask_png,
                max_threads=args.max_threads or config.agent_tools.max_threads,
                threshold=args.threshold,
                keep_mode=args.keep_mode,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        elif args.command == "recommend-preprocess":
            from .preprocessing_recommender import recommend_preprocessing_plan

            result = recommend_preprocessing_plan(
                args.input,
                biological_metadata_json=args.biological_metadata_json,
                biological_metadata_path=args.biological_metadata_path,
                max_threads=args.max_threads or config.agent_tools.max_threads,
                sample_size=args.sample_size or config.agent_tools.qc_sample_size,
                noise_method=args.noise_method or config.agent_tools.noise_method,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        elif args.command == "doctor":
            from .doctor import run_self_check

            result = run_self_check(
                args.input,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
                config_path=args.config,
                require_mcp=args.require_mcp,
            )
        elif args.command == "usage-modes":
            usage = build_usage_modes_payload(config)
            if args.json:
                sys.stdout.write(json.dumps(usage, indent=2, ensure_ascii=False))
                sys.stdout.write("\n")
            else:
                sys.stdout.write(render_usage_modes_text(usage))
                sys.stdout.write("\n")
            return 0
        elif args.command == "validate":
            from .validation_tools import validate_integrated_conclusion

            operations = json.loads(args.operations_json) if args.operations_json else None
            result = validate_integrated_conclusion(
                args.input,
                run_preprocessing_first=not args.no_preprocess,
                preprocessing_operations=operations,
                biological_metadata_json=args.biological_metadata_json,
                biological_metadata_path=args.biological_metadata_path,
                chemical_library_dir=args.chemical_library_dir,
                chemical_source_name=args.chemical_source_name,
                chemical_ion_mode=args.chemical_ion_mode,
                chemical_ppm_tolerance=args.chemical_ppm_tolerance,
                chemical_max_mz_values=args.chemical_max_mz_values,
                chemical_max_candidates_per_mz=args.chemical_max_candidates_per_mz,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
        else:
            from .workflow_tools import analyze_imzml_end_to_end

            operations = json.loads(args.operations_json) if args.operations_json else None
            result = analyze_imzml_end_to_end(
                args.input,
                run_preprocessing_first=not args.no_preprocess,
                preprocessing_operations=operations,
                biological_metadata_json=args.biological_metadata_json,
                biological_metadata_path=args.biological_metadata_path,
                chemical_library_dir=args.chemical_library_dir,
                chemical_source_name=args.chemical_source_name,
                chemical_ion_mode=args.chemical_ion_mode,
                chemical_ppm_tolerance=args.chemical_ppm_tolerance,
                chemical_max_mz_values=args.chemical_max_mz_values,
                chemical_max_candidates_per_mz=args.chemical_max_candidates_per_mz,
                output_dir=args.output_dir or config.agent_tools.default_output_dir,
            )
    except ModuleNotFoundError as exc:
        sys.stderr.write(
            "Missing runtime dependency while loading the selected MassFlow tool: "
            f"{exc.name}\n"
        )
        return 2

    sys.stdout.write(json.dumps(make_json_safe(result.to_dict()), indent=2, ensure_ascii=False))
    sys.stdout.write("\n")
    return 0 if result.ok else 1


def build_usage_modes_payload(config) -> dict[str, object]:
    project_dir = str(Path.cwd())
    return {
        "formal_mode": {
            "name": "Claude Code automatic MCP mode",
            "description": "Use Claude Code with stdio MCP config. Do not start the server manually.",
            "startup": [
                f'cd "{project_dir}"',
                "claude",
            ],
            "mcp_transport": config.mcp_server.transport,
        },
        "manual_debug_mode": {
            "name": "Manual local debug mode",
            "description": "Start the MCP server yourself only for debugging HTTP transport.",
            "startup": [
                f'cd "{project_dir}"',
                (
                    "massflow-mcp-server --transport streamable-http "
                    f"--host {config.mcp_server.host} --port {config.mcp_server.port}"
                ),
            ],
            "mcp_transport": "streamable-http",
        },
        "switch_rule": "Formal use means Claude Code auto-starts MCP. Manual debug means you start the HTTP server yourself.",
    }


def _add_biological_metadata_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--biological-metadata-json",
        help="Structured biological metadata as a JSON object string.",
    )
    parser.add_argument(
        "--biological-metadata-path",
        help="Path to a JSON file containing structured biological metadata.",
    )


def _add_chemical_annotation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--chemical-library-dir",
        help="Directory containing local chemical library TSV/CSV files.",
    )
    parser.add_argument(
        "--chemical-source-name",
        default="auto",
        help="Source label used in chemical annotation results.",
    )
    parser.add_argument(
        "--chemical-ion-mode",
        choices=["positive", "negative"],
        help="Override ion mode for local chemical annotation.",
    )
    parser.add_argument(
        "--chemical-ppm-tolerance",
        type=float,
        default=5.0,
        help="PPM tolerance used for local chemical annotation.",
    )
    parser.add_argument(
        "--chemical-max-mz-values",
        type=int,
        default=12,
        help="Maximum representative m/z values extracted from QC for annotation.",
    )
    parser.add_argument(
        "--chemical-max-candidates-per-mz",
        type=int,
        default=5,
        help="Maximum candidate annotations returned for each m/z value.",
    )


def render_usage_modes_text(payload: dict[str, object]) -> str:
    formal_mode = payload["formal_mode"]
    manual_mode = payload["manual_debug_mode"]
    return (
        "Formal use:\n"
        f"- {formal_mode['description']}\n"
        f"- Command 1: {formal_mode['startup'][0]}\n"
        f"- Command 2: {formal_mode['startup'][1]}\n"
        "Manual debug:\n"
        f"- {manual_mode['description']}\n"
        f"- Command 1: {manual_mode['startup'][0]}\n"
        f"- Command 2: {manual_mode['startup'][1]}\n"
        f"Switch rule: {payload['switch_rule']}"
    )


if __name__ == "__main__":
    raise SystemExit(main())
