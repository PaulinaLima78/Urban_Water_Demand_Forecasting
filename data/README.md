# Data

Four input files are required. The script locates them by pattern, so exact
file names may vary provided the identifying token is present.

| Pattern matched | Expected content |
|---|---|
| `agua*quito*`, `demand*` | Daily treated flow and calendar flags |
| `*bellavista*precip*`, `c05*` | Daily rainfall, Bellavista gauge |
| `*inaquito*precip*`, `p09*` | Daily rainfall, Inaquito gauge |
| `oni*` | Monthly Oceanic Nino Index |

File names use ASCII characters only, to avoid encoding problems across
operating systems.

## Columns

**Demand** — `n`, `date` (d/m/Y), `day-week`, `day-month`, `val` (L/s),
`weekend`, `holiday`, `sequia`.

Note on `sequia`: this column is a deterministic June–September indicator,
identical in every year. It marks the dry SEASON and is not a drought index;
it is collinear with the month factor and its main effect is therefore not
identifiable. The script reports this explicitly.

**Precipitation** — `fecha` (Y/m/d), `valor` (mm), completeness flags.

**ONI** — `YR`, `MON`, `TOTAL C`, `limAdjust`, `ANOM` (the anomaly is used).
