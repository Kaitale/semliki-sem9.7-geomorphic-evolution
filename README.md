# Meander migration, cut-off risk and bed stability on the lower Semliki River, Albertine Rift: 26 years of Landsat and Sentinel-2 and a 75-year gauging record

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24367442.svg)](https://doi.org/10.5281/zenodo.24367442) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Kaitale/semliki-sem9.7-geomorphic-evolution/blob/main/reproduce_paper.ipynb)

**Monica Nafula (KAI Water & Energy); Charles Kakooza (Ministry of Water and Environment)**

Preprint (not peer reviewed), version 2.0, 19 November 2026, prepared to the author guidelines of *Journal of Hydrology: Regional Studies*.

## Abstract

**Study region:** The lower Semliki River on the floor of the Albertine Rift, from 10 km above the Bweramule gauge (Uganda–Democratic Republic of Congo border) to Lake Albert, about 102 km of a 9,917 km² basin draining Lake Edward and the Rwenzori Mountains. The reach floods Ntoroko District, Uganda, from the river and from Lake Albert backwater.

**Study focus:** Flood models of the reach treat the channel as fixed. We measure how it has moved since 2000 and whether its bed has changed, using 27 annual Landsat and 11 Sentinel-2 low-flow water masks (MNDWI with an Otsu threshold), linear-regression bank rates with 95 % confidence intervals on 4,090 transects at 25 m spacing, bend-scale statistics and 50-year setbacks for 125 bends, a search for meander necks, a specific-gauge analysis of 191 discharge gaugings (1950–2025) and low-water trends at satellite-altimetry virtual stations.

**New hydrological insights for the region:** Between 2000–04 and 2020–26, 265 ha of floodplain became channel and 216 ha of channel became land in the river reach; two-fifths of today’s low-flow channel lies on land that was dry 20 years earlier. Outer banks of the typical bend retreat at 4.5–5.0 m/yr in the meander belt, 1.5 m/yr between km 40 and 60 and about 1 m/yr in the lake-controlled reach; seven bends exceed 10 m/yr. Scaled by width (0.056 widths/yr), migration is three to five times faster than reported for vegetated rivers worldwide. A 97 m neck at km 14.2 is closing at about 24 m/yr and is likely to cut off within 3–10 years. Stage for a given flow at the gauge has no net trend since 1950, although the control rose about 0.35 m in 1955–1975 and fell back; rising low water in the lower reach tracks the 2.0 m rise of Lake Albert, not sediment. Median 50-year setbacks are 205 m, rising to 578 m for nine bends in ten.

## Reproduce the paper

One notebook, top to bottom:

```bash
pip install -r requirements.txt
jupyter lab reproduce_paper.ipynb      # then Run All
```

or open it in [Google Colab](https://colab.research.google.com/github/Kaitale/semliki-sem9.7-geomorphic-evolution/blob/main/reproduce_paper.ipynb); the notebook downloads anything it needs from this repository.

## Contents

| File | What it is |
|---|---|
| `Nafula_Kakooza_2026_Semliki_meander_migration_preprint.pdf` | the paper |
| `reproduce_paper.ipynb` | the one notebook: recomputes Tables 2–5 and Figs. 2, 4–7 of the paper (about 2 min) |
| `semliki_bankmig.py` | the transect, bank-rate, bend and neck analysis (imported by the notebook) |
| `data_masks_and_results.zip` | 27 Landsat and 11 Sentinel-2 annual low-flow water masks (.npz) and the per-transect, bend and neck results of the three records (.csv) |
| `extras/build_water_masks_gee.ipynb` | optional: rebuilds the annual masks in Google Earth Engine (needs a GEE account) |
| `requirements.txt` | numpy, scipy, matplotlib (plus earthengine-api for the extra notebook) |

## Data not included

The Bweramule discharge gaugings and the HEC-RAS model are held by the Ministry of Water and Environment, Uganda, and are available on request. Landsat imagery courtesy of the U.S. Geological Survey; contains modified Copernicus Sentinel data (2016–2026); altimetry from Hydroweb (CNES/Theia) and DAHITI (DGFI-TUM).

## How to cite

Nafula, M., Kakooza, C., 2026. Meander migration, cut-off risk and bed stability on the lower Semliki River, Albertine Rift: 26 years of Landsat and Sentinel-2 and a 75-year gauging record. Preprint, version 2.0. Zenodo. https://doi.org/10.5281/zenodo.24367442

## Licence

Manuscript and data: CC BY 4.0. Code: MIT. See `LICENSE`.

## Funding

Carried out under the Semliki Flood Forecasting System of the Ministry of Water and Environment, Uganda, with support from the World Food Programme.

## Version history

* **2.0** (19 Nov 2026): the work as a journal-style preprint with one reproduction notebook.
* 1.x (Oct 2026): the original technical report and analysis notebooks (superseded).
