# MPTRAC Web Docker Deployment

This directory contains the Docker Compose configuration for deploying the MPTRAC web application using the prebuilt Docker image.

## Requirements

* Docker and Docker Compose
* The meteorological data filesystem mounted on the host

## Configuration

Create the local environment file:

```bash
cp .env.example .env
```

Edit `.env` as needed:

```bash
MPTRAC_VERSION=0.1
MPTRAC_PORT=80
DATA_HOST_PATH=/mnt/slmetdata_mnt
MET_ACCESS_PATH=/mnt/slmetdata_mnt/met_data/
```

`DATA_HOST_PATH` must point to the host directory containing the mounted meteorological data.

## Deployment

Start the service:

```bash
docker compose pull
docker compose up -d
```

Check the container:

```bash
docker compose ps
docker compose logs
```

Stop the service:

```bash
docker compose down
```

The web application is available on the configured `MPTRAC_PORT`, for example:

```text
http://<server-address>:80/
```
