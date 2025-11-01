# Building Docker Image for x86_64 on ARM Mac

## First-time Setup

```bash
docker buildx create --name multiplatform --use
docker buildx inspect --bootstrap
```

## Build x86_64 Image

```bash
docker buildx build --platform linux/amd64 -t phlawd:latest .
```

To build and load into your local Docker (for testing on your Mac):

```bash
docker buildx build --platform linux/amd64 -t phlawd:latest --load .
```

To build and push to a registry:

```bash
docker buildx build --platform linux/amd64 -t phlawd:latest --push .
```
