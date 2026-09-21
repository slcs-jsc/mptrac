# MPTRAC Web Runner — Docker HTTPS Deployment

This directory provides an HTTPS deployment of the MPTRAC Web Runner using Docker Compose and Nginx.

Nginx acts as the public HTTPS reverse proxy, while the MPTRAC Web Runner communicates with Nginx over the internal Docker network using HTTP.

```text
Client
  |
  | HTTPS
  v
Nginx :443
  |
  | HTTP (Docker network)
  v
MPTRAC Web Runner :80
```

## Requirements

Docker with the Docker Compose plugin must be installed on the host.

A TLS certificate and corresponding private key must already be available on the host. Certificate acquisition and renewal are intentionally outside the scope of this deployment.

The deployment is certificate-authority independent. Certificates issued by Let's Encrypt, HARICA, an institutional CA, or another trusted CA can be used.

The MPTRAC meteorological data must also be available on the Docker host.

## Configuration

Create the local environment file from the provided example:

```bash
cp .env.example .env
```

Then edit `.env` for the target system.

Example:

```dotenv
MPTRAC_WEB_VERSION=0.1
MPTRAC_DOMAIN=example.org
MPTRAC_PORT=443

DATA_HOST_PATH=/path/to/data
MET_ACCESS_PATH=/path/to/data/met_data/

TLS_FULLCHAIN_CERT_PATH=/path/to/fullchain.pem
TLS_PRIVATE_KEY_PATH=/path/to/privatekey.pem
```

`MPTRAC_WEB_VERSION` specifies the MPTRAC Web Runner Docker image version.

`MPTRAC_DOMAIN` is the hostname served by Nginx and must correspond to the TLS certificate.

`MPTRAC_PORT` specifies the public HTTPS port. Port `443` is recommended for a standard HTTPS deployment.

`DATA_HOST_PATH` specifies the data directory on the Docker host that is mounted read-only into the MPTRAC container.

`MET_ACCESS_PATH` specifies the meteorological-data location used by the MPTRAC Web Runner.

`TLS_FULLCHAIN_CERT_PATH` specifies the host path to the server certificate including the required certificate chain.

`TLS_PRIVATE_KEY_PATH` specifies the corresponding TLS private key.

Certificate and private-key files are mounted read-only into the Nginx container and are not stored in the MPTRAC repository.

## Validate the configuration

Before starting the deployment, inspect the resolved Docker Compose configuration:

```bash
docker compose config
```

Required variables and required bind-mount paths are validated by the Compose configuration.

## Start the service

```bash
docker compose up -d
```

Check the running containers:

```bash
docker compose ps
```

Validate the generated Nginx configuration:

```bash
docker exec nginx-mptrac-https nginx -t
```

The service should then be available at:

```text
https://<MPTRAC_DOMAIN>/
```

If a non-standard `MPTRAC_PORT` is configured, include it explicitly:

```text
https://<MPTRAC_DOMAIN>:<MPTRAC_PORT>/
```

## Stop the service

```bash
docker compose down
```

## TLS certificate management

This deployment consumes existing TLS certificate files but does not obtain or renew certificates.

Certificate issuance and renewal should be handled independently using the mechanism appropriate for the deployment environment, such as Let's Encrypt/ACME, HARICA, or an institutional certificate-management service.

After replacing or renewing certificate files, Nginx may need to be reloaded or restarted so that it uses the updated certificate.
