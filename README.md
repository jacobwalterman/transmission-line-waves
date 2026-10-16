# transmission-line-waves

A nanobind-wrapped C++ solver for a 4-segment transmission-line wave equation, with Python scripts for plotting.

![Example output](docs/example.png)

The command-line interface is being reworked, so run instructions are coming soon. For now, this page covers setup only.

## Install with pip

Needs Python 3.12+, a C++ compiler, and CMake. From a checkout of this repo:

```
pip install .
```

## Docker

Build the image:

```
docker build -t transmission-line-waves .
```

## Local development (Nix)

Requires [Nix](https://nixos.org/) and [direnv](https://direnv.net/).

```
direnv allow
```

This enters the flake's dev shell (compiler, CMake, Eigen, Python) and creates and activates a project-local venv automatically. Then, once per checkout and again after changing any C++ source:

```
pip install .
```
