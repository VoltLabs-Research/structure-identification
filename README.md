# StructureIdentification

`StructureIdentification` is the shared reconstructed-state and cluster-graph library used by CNA, PTM, PSM, ElasticStrain and OpenDXA.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- StructureIdentification
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Output

This plugin builds a static library (`libstructure-identification.a`) under `build/build/Release/`. It does not produce a standalone CLI; downstream plugins (CNA, PTM, OpenDXA, etc.) consume it through Conan.

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

## CLI

`StructureIdentification` does not install a standalone CLI executable.
