# Optional workflows

Install the project and activate its Python environment before running these scripts. They are not needed to install or use the CLI.

- `generate_featured_molecules*.sh`: illustrative PubChem molecule batches.
- `validate_sugar_batch.py`, `validate_flatness_models.py`: larger geometry validation runs that write to `out/`.
- `prepare_*`, `generate_*`, `index_*`: collection-specific workflows. Many expect manifests or previous models in `out/`; those local collections are intentionally not distributed.
- `audit_*`, `repair_*`, `refresh_*`, `organize_*`, `rerender_*`, `report_*`, `check_*`, and collection validators: maintenance tools for existing collections. Read each script before running; some move or replace existing outputs.

For a self-contained, offline starting point, use `3dprintase batch examples/sugars.json --no-render --output out/example` instead. The CLI/API is the supported general interface; these scripts record individual workflows and may require local inputs.
