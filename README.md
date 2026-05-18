# gdml-to-mesh

A tool for converting GDML (Geometry Description Markup Language) files into mesh formats for simulation and visualization purposes, with advanced support for boundary surface extraction using Open CASCADE (OCC).

## Features

- Parses GDML geometry files.
- Extracts and processes geometrical shapes and **boundary surfaces**.
- Uses Open CASCADE (pythonOCC) for CAD processing and mesh extraction.
- Generates mesh files compatible with physics simulations and 3D visualization.
- Supports custom configuration for mesh and geometry extraction.
- Enforces code quality with automated `pre-commit` hooks.

## Prerequisites

1. **Python 3.x**
   Install Python 3 from [python.org](https://www.python.org/downloads/).

2. **Python Virtual Environment (Optional Recommended)**
   Manage dependencies in an isolated environment.

3. **Open CASCADE (OCC) / pythonOCC**
   Used for CAD geometry processing and boundary surface extraction.
