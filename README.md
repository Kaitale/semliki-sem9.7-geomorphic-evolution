# Geomorphic Evolution of the Lower Semliki: Bed and Bank Change, 2000–2026

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24367443.svg)](https://doi.org/10.5281/zenodo.24367443)

**SFFS Technical Report SEM.9.7, version 1.1 (4 October 2026)** — Semliki Flood Forecasting System, Ministry of Water and Environment, Uganda.

The lower Semliki — the 92 km of river between the Bweramule gauge and Lake Albert — is the reach that floods Ntoroko District. This report measures how its banks and bed have changed since 2000 from 27 annual Landsat and 11 Sentinel-2 low-flow water masks, 4,090 transects, 139 bends, the Bweramule gaugings (1950–2025) and satellite altimetry, and says what the findings mean for the flood forecast and for planning.

**Key findings:** outer banks retreat at 4.5–5.0 m/yr at the typical bend of the first 40 km below the gauge (the fastest bends 11–13 m/yr) and about 1 m/yr beyond km 60; a 97 m meander neck at km 14.2 is the next cut-off (roughly 3–10 years at present rates); the bed at the gauge shows no net aggradation in 70 years; rising low water in the lower reach is Lake Albert, not sediment; a 50-year setback of about 205 m keeps half of the bends clear and about 580 m keeps nine in ten clear.

## Files

| File | What it is |
|---|---|
| `SEM.9.7.R.1.1.SFFS.Geomorphic.Analysis.Bed.and.Bank.Evolution.pdf` | The technical report |
| `SEM.9.7_data_masks_and_results.zip` | `data/masks/`: 27 Landsat (2000–2026) and 11 Sentinel-2 (2016–2026) annual low-flow channel-water masks, `.npz` (`water`, `valid`, `transform`), UTM zone 35N (EPSG:32635). `data/results/`: the published transect, bend, neck and bank-position tables of the three records, and the channel area by year |
| `semliki_bankmig.py` | Transects, bank positions, retreat rates with 95 % intervals, bends, necks and setbacks |
| `01_water_masks_gee.ipynb` | Builds the annual masks in Google Earth Engine (needs an Earth Engine account; shown un-run) |
| `02_bank_retreat_bends_necks.ipynb` | Re-runs the transect analysis for the three records and checks it against the published results |
| `03_planform_change_and_width.ipynb` | Epoch composites, erosion and deposition by zone, channel width by year |
| `04_bed_specific_gauge_altimetry_methods.ipynb` | Specific-gauge and altimetric low-water methods, run on a labelled synthetic example |

Run with Python 3.10+ and `pip install -r requirements.txt`; open the notebooks from this folder. The zip is extracted automatically on first run.

**Data not included.** The discharge gaugings and gauge levels at Bweramule (DWRM 85205) and the HEC-RAS model are held by the Ministry of Water and Environment, Uganda, and are available from the Ministry on request. Altimetric water levels are available from Hydroweb (Theia/CNES). Contains modified Copernicus Sentinel data (2016–2026); Landsat imagery courtesy of the U.S. Geological Survey.

## Authors

Monica Nafula and Charles Kakooza.

Prepared under the Semliki Flood Forecasting System (SFFS) of the Ministry of Water and Environment (MWE), Uganda, with support from the World Food Programme (WFP).

## Licence

Report and data: CC BY 4.0. Code and notebooks: MIT. See `LICENSE`.

## How to cite

Zenodo record: https://doi.org/10.5281/zenodo.24367443 (report, data and notebooks). See also `CITATION.cff`.
