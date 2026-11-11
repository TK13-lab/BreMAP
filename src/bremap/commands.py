"""Add one-run orchestration without changing the frozen historical dispatcher."""
import argparse
import sys
from zipfile import BadZipFile

from . import cli
from .workflow import run_workflow


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] == 'run':
        parser = argparse.ArgumentParser(prog='bremap run', description='Execute a reproducible JSON workflow recipe')
        parser.add_argument('--recipe', required=True)
        parser.add_argument('--out', help='Override recipe output directory')
        parser.add_argument('--overwrite', action='store_true', help='Replace recognized outputs only after success')
        options = parser.parse_args(args[1:])
        try:
            return run_workflow(options.recipe, out_dir=options.out, overwrite=options.overwrite)
        except (ValueError, OSError, KeyError, TypeError, BadZipFile) as exc:
            parser.error(str(exc))
    if not args or args[0] in ('-h', '--help'):
        print('One-run workflow: bremap run --recipe recipe.json [--out DIRECTORY] [--overwrite]\n')
    return cli.main(args)


def console_main():
    main()
