# Process  barcoded sequencing libraris with UMIs

This script processes libraries that contain cell barcodes and UMIs in one read, and 3' sequences in the other. This script parses read1, generating cell-specific fastq files for read2 with UMIs appended to the read headers in read2. These Each demultiplexed sample in then aligned to a reference genome using STAR (Dobin *et al.*, 2013) and PCR duplicates are removed using UMI-tools (Smith *et al.* 2017).

## Getting Started

For usage, use

`process_800_cell_chip.py -h`

### Prerequisites

* python (developed on python 2.7.8)
* UMI-tools ()
* STAR (developed using STAR 2.5.0b)
* pandas
* HTSeq-count
* samtools

Additionally, an indexed reference genome in STAR format is required as well as the GTF that was used to create it.

### Installing

Currently, just run the script and ensure the dependencies above are satisfied.

## Running the tests

Limited testing has been performed, but it requires files that are too large for github. This will be remedied later.

## To do:

* Create tests with small memory fingerprint
* fill out function documentation
* Extend logging to keep track of read counts

## Authors

* **Rob Carter** - (r.andrew.carter@gmail.com)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to - - for the mapping between cell barcode 6mer and row number

## References
	Smith T, Heger A, Sudbery I. 2017. UMI-tools: modeling sequencing errors in Unique
	Molecular Identifiers to improve quantification accuracy. Genome Res. 27:491-499
	
	Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, 
	Gingeras TR. 2013. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 29:15-21
	
