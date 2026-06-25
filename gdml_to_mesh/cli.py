"""
cli.py — gdml-to-mesh command line interface.

Usage:
    gdml-to-mesh run detector.gdml
    gdml-to-mesh run detector.gdml --output-dir ./output
    gdml-to-mesh build
    gdml-to-mesh check
    gdml-to-mesh visualize
    gdml-to-mesh visualize --surface detector
"""

import argparse
import sys
from pathlib import Path


def cmd_run(args):
    from gdml_to_mesh import run
    result = run(
        args.gdml,
        output_dir=args.output_dir,
        rebuild=args.rebuild,
        verbose=args.verbose,
        jobs=args.jobs,
    )
    print(result.summary())


def cmd_build(args):
    from gdml_to_mesh import _build
    _build.build(verbose=args.verbose, jobs=args.jobs)


def cmd_check(args):
    from gdml_to_mesh import check_dependencies, _build
    check_dependencies()
    print(f"\n  binary built : {_build.is_built()}")
    if _build.is_built():
        print(f"  binary path  : {_build.binary_path()}")


def cmd_visualize(args):
    import subprocess
    viz = Path(__file__).parent.parent / "geant4_extract" / "visualize.py"
    if not viz.exists():
        print(f"Visualizer not found at {viz}")
        sys.exit(1)
    cmd = [sys.executable, str(viz)]
    if args.surface:
        cmd += ["--surface", args.surface]
    if args.no_volumes:
        cmd += ["--no-volumes"]
    if args.list:
        cmd += ["--list"]
    subprocess.run(cmd)


def main():
    parser = argparse.ArgumentParser(
        prog="gdml-to-mesh",
        description="LEGEND-Theia GDML geometry engine"
    )
    sub = parser.add_subparsers(dest="command")

    # run
    p_run = sub.add_parser("run", help="Run the geometry engine on a GDML file")
    p_run.add_argument("gdml", help="Path to .gdml file")
    p_run.add_argument("--output-dir", default=None, help="Output directory (default: cwd)")
    p_run.add_argument("--rebuild", action="store_true", help="Force rebuild of C++ binary")
    p_run.add_argument("--verbose", action="store_true")
    p_run.add_argument("--jobs", type=int, default=8)

    # build
    p_build = sub.add_parser("build", help="Build the C++ binary")
    p_build.add_argument("--verbose", action="store_true")
    p_build.add_argument("--jobs", type=int, default=8)

    # check
    sub.add_parser("check", help="Check dependency paths and build status")

    # visualize
    p_viz = sub.add_parser("visualize", help="Launch the geometry visualizer")
    p_viz.add_argument("--surface", choices=["blackbody", "specular", "detector", "default"], default=None)
    p_viz.add_argument("--no-volumes", action="store_true")
    p_viz.add_argument("--list", action="store_true")

    args = parser.parse_args()

    if args.command == "run":
        cmd_run(args)
    elif args.command == "build":
        cmd_build(args)
    elif args.command == "check":
        cmd_check(args)
    elif args.command == "visualize":
        cmd_visualize(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
