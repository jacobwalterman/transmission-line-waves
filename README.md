# transmission-line-waves

A nanobind-wrapped C++ solver for a 4-segment transmission-line wave equation, with Python driver scripts for plotting.

## Local development (Nix)

Requires [Nix](https://nixos.org/) and [direnv](https://direnv.net/).

```bash
direnv allow
```

This enters the flake's dev shell (compiler, CMake, Eigen, Python) and auto-creates/activates a project-local venv via direnv's layout python3 — no manual `source .venv/bin/activate` needed. Then, once per checkout (and again after changing any C++ source):

```bash
pip install .
```

Run the solver:

```bash
python python/run_transmission_line.py [box|flat|sin] [out_path.png]
# defaults to "box" and transmission_line.png
```

`python/run_engine.py` is a separate, interactive (`plt.show()`) demo script for the simpler single-segment solver — local-only, not part of the Docker image.

## Docker

Build:

```bash
docker build -t transmission-line-waves:dev .
```

Run and retrieve the output PNG via `docker cp` (avoids bind-mount permission issues, e.g. SELinux on Fedora):

```bash
sudo docker run --name tlw-run transmission-line-waves:dev
sudo docker cp tlw-run:/app/transmission_line.png ./transmission_line.png
sudo docker rm tlw-run
```

Override the initial condition / output filename via CMD args:

```bash
sudo docker run --name tlw-run transmission-line-waves:dev sin /app/custom.png
sudo docker cp tlw-run:/app/custom.png ./custom.png
sudo docker rm tlw-run
```

If you'd rather use a live bind mount instead of `docker cp`, add `:z` on SELinux-enforcing systems (e.g. Fedora) so the container can actually write to it:

```bash
mkdir -p /tmp/out
sudo docker run --rm -v /tmp/out:/app/out:z transmission-line-waves:dev box /app/out/box.png
```
