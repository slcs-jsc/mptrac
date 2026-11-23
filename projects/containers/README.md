# MPTRAC containers

This directory contains Docker-based deployment configurations for the
[MPTRAC Web Runner](../web_runner/). The Docker setup uses the image
`amirnikfal/mptrac-web` provided by Amir Nikfal and supports HTTP and
HTTPS deployments.

## Directory structure

- [`docker/Dockerfile.web`](docker/Dockerfile.web) defines the Web Runner image.
- [`docker/Dockerfile.web.dockerignore`](docker/Dockerfile.web.dockerignore)
  limits the Docker build context to the required files.
- [`docker/http/`](docker/http/) contains a simple HTTP deployment using
  Docker Compose.
- [`docker/https/`](docker/https/) contains an HTTPS deployment with
  Nginx as a reverse proxy and externally managed TLS certificates.

Both deployments use the prebuilt Docker image
`amirnikfal/mptrac-web:0.1` by default; their image tags can be changed
through their respective environment files.

## Building the image

The Dockerfile expects the MPTRAC repository root as its build context.
The MPTRAC executables `src/atm_init`, `src/trac`, and `src/atm_conv`
must already be available: the Dockerfile copies these binaries rather
than compiling MPTRAC. They must be compatible with Linux/amd64 and the
container runtime.

From the repository root, run:

```bash
docker build -f projects/containers/docker/Dockerfile.web \
    -t mptrac-web:local .
```

The image is based on Python 3.12 (Debian Bookworm), installs the Web
Runner's Python dependencies and Gunicorn, copies the required MPTRAC
binaries and model data, and preloads Cartopy map data. It runs as a
non-root user and starts Gunicorn on container port 80.

To use this local image with Docker Compose, adjust the `image:`
reference in the chosen Compose file accordingly. The default Compose
configurations instead pull Amir's published image.

## Deployment

Choose one of the configurations below and follow its README for the
required environment variables, data mounts, and startup commands:

- **[HTTP deployment](docker/http/README.md):** Directly publishes the
  Web Runner's HTTP port. Useful for local or controlled environments.
- **[HTTPS deployment](docker/https/README.md):** Places Nginx in front of
  the Web Runner, terminates TLS, and exposes an HTTPS port. This setup
  requires a hostname and existing certificate/key files.

For either variant, meteorological data must be accessible on the
Docker host and mounted read-only into the container. Configure
`DATA_HOST_PATH` and `MET_ACCESS_PATH` as described in the deployment
README. The HTTPS setup additionally needs TLS certificate paths.

For example, to start the HTTP variant:

```bash
cd projects/containers/docker/http
cp .env.example .env
# Edit .env to match the host and data paths.
docker compose up -d
docker compose ps
```

Use `docker compose logs` to inspect the service and
`docker compose down` to stop it.

For application details outside container deployment, see the
[MPTRAC Web Runner](../web_runner/).
