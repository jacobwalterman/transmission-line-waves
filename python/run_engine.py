import subprocess
from pathlib import Path

EXE_PATH = Path(__file__).resolve().parent.parent / "build" / "coupled_piezo_heat"
def run_engine(n=10, freq=1.0, exe_path=EXE_PATH):
    result = subprocess.run(
        [exe_path, "--n", str(n), "--freq", str(freq)],
        capture_output=True,
        text=True,
        check=True,
    )
    values = [float(line) for line in result.stdout.strip().splitlines()]
    return values

if __name__ == "__main__":
    values = run_engine(n=20, freq=2.0)
    print(values)
