# StructureIdentification

`StructureIdentification` is the shared reconstructed-state and cluster-graph library used by CNA, PTM, PSM, ElasticStrain and OpenDXA.

## CLI

`StructureIdentification` does not install a standalone CLI executable.

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/StructureIdentification
conan create . -nr
```
