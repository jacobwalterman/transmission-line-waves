# transmission-line-waves

A solver for a transmission-line model made of four wave equations joined end to end. The junctions between segments can be continuous or discontinuous, and controllers act at the discontinuities. This is the numerical implementation of the model from [our paper](https://www.aimsciences.org/article/doi/10.3934/eect.2026007) in *Evolution Equations and Control Theory* (2026), which gives the full model and its application to transmission lines.

The numerical core is written in C++ and wrapped with nanobind. Python handles configuration and plotting through a small command-line tool.

![Example output](docs/example.png)

*Left to right: positions, velocities, and normalized energy E(t)/E(0) over time. The t = 0 frame is the exact initial condition rather than the computed solution, which is why the energy jumps at the start.*

## Quick start (Docker)

You only need Docker. Build the image once:

```
docker build -t transmission-line-waves .
```

Then run it in any folder:

```
docker run --rm -v "$PWD:/out:z" transmission-line-waves box --out box.png
```

This saves `box.png` in your current folder. The first argument is the initial condition (`box`, `flat`, or `sin`, default `box`). The arguments are the same as the `tlw` command below. The `-v` part is what makes the output show up on your machine, so don't leave it out.

The `:z` is for SELinux systems like Fedora and is harmless elsewhere. Add `sudo` if your Docker setup needs it.

## Install with pip

Needs Python 3.12+, a C++ compiler, and CMake.

```
pip install git+https://github.com/jacobwalterman/transmission-line-waves
```

or, from a checkout of this repo:

```
pip install .
```

Then:

```
tlw box --out box.png
tlw --help
```

## Options

```
tlw [-h] [--out OUT] [--config CONFIG] [{box,flat,sin}]
```

- `--out`: output PNG path (default `transmission_line.png`)
- `--config`: path to your own params TOML file (default: the built-in one)

## Local development (Nix)

Requires [Nix](https://nixos.org/) and [direnv](https://direnv.net/).

```
direnv allow
```

This enters the flake's dev shell (compiler, CMake, Eigen, Python) and creates and activates a project-local venv automatically. Then, once per checkout and again after changing any C++ source:

```
pip install .
```

## License

MIT. See [LICENSE](LICENSE).

## Troubleshooting

If the bind mount gives you permission errors, copy the output out of the container instead:

```
docker run --name tlw-run transmission-line-waves sin --out custom.png
docker cp tlw-run:/out/custom.png ./custom.png
docker rm tlw-run
```
