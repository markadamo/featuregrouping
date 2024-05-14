# featuregrouping
Group MS1 features for mTRAQ labeled peptides

Install dependencies with:
`pip3 install -r requirements.txt`

Example command to test on provided data (write to stdout):
`python3 features.py data/ffc.featureXML -p data/psms.tsv`

```
usage: python3 features.py [-h] [-t {featureXML,biosaur2}] [-p PSM_FILE] [-o OUTPUT_FILE]
                   [-ppm MZ_TOL_PPM] [-rt RT_TOL] [-nt MAX_TAGS_PER_PEPTIDE]
                   feature_file

mTRAQ feature grouping

positional arguments:
  feature_file

options:
  -h, --help            show this help message and exit
  -t {featureXML,biosaur2}, --input_type {featureXML,biosaur2}
                        Specify 'biosaur2' if the input file is a biosaur2 tsv
  -p PSM_FILE, --psm_file PSM_FILE
                        Optionally provide a PSM tsv file to perform QC. Adds matched
                        peptide sequences to output.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Output path for tsv file. Omit to output to stdout.
  -ppm MZ_TOL_PPM, --mz_tol_ppm MZ_TOL_PPM
                        M/z tolerance in ppm.
  -rt RT_TOL, --rt_tol RT_TOL
                        Retention time tolerance in seconds.
  -nt MAX_TAGS_PER_PEPTIDE, --max_tags_per_peptide MAX_TAGS_PER_PEPTIDE
                        Maximum number of mTRAQ tags to consider per peptide.
```
