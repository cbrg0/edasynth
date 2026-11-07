"""This file lets the program be started with the command python -m edasynth. It only
passes control to the file cli.py, which knows all the commands.
"""

from .cli import main

raise SystemExit(main())
