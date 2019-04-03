README.vhs lists some statistics; e.g. number of tiles, framesets, sources

```
usage: vhsph3_getdata_vsa.py [-h] [-t] [--debug] [-l] [-i INPUT] [-q QUERY]
                             [-f FORMAT] [-o OUTPUT] [--config CONFIG]
                             [--outpath OUTPATH] [-e EMAIL] [-d DB]
                             [--dr RELEASE]

Download VHS data tile by tile from VSA

optional arguments:
  -h, --help            show this help message and exit
  -t, --test            run test (default: False)
  --debug               run in debug mode (default: False)
  -l, --limit           limit SQL to 100 rows with TOP (default: False)
  -i INPUT, --input INPUT
                        input metadata file (default: None)
  -q QUERY, --query QUERY
                        specify query on the command line (default: None)
  -f FORMAT, --format FORMAT
                        Output format (csv/fits/fits.gz/html) (default:
                        .fits.gz)
  -o OUTPUT, --output OUTPUT
                        Output filename (default: None)
  --config CONFIG       Specify the config file (default: ./)
  --outpath OUTPATH     Output path (default: ./)
  -e EMAIL, --email EMAIL
                        Email address for results of long queries; should not
                        be used when batch processing is required (default: )
  -d DB, --db DB        VSA Database to use (default: VHSv20160507)
  --dr RELEASE          VSA Data release to use (default: VHSDR2)

Work in progress
```
