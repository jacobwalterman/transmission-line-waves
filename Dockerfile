# ---- Builder: compile the extension and build a wheelhouse ----
FROM python:3.12-slim AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY pyproject.toml CMakeLists.txt ./
COPY src ./src

RUN pip install --no-cache-dir --upgrade pip \
    && pip wheel --no-cache-dir -w /wheels .

# ---- Runtime: slim image with only the installed package + script ----
FROM python:3.12-slim AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir --no-index --find-links=/wheels transmission_line_waves \
    && rm -rf /wheels

COPY python/run_transmission_line.py ./run_transmission_line.py

ENTRYPOINT ["python", "run_transmission_line.py"]
CMD ["box", "transmission_line.png"]
