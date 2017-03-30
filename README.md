# Process the 800 cell IFC from Fluidigm

The fluidigm 800-cell IFC for the C1 has both barcodes and UMIs in read 1, whereas read 2 contains the actual sequence of a polyadenylated RNA. This script parses read1, generating cell-specific fastq files for read2 with UMIs appended to the read headers in read2. These Each demultiplexed sample in then aligned to a reference genome using STAR () and PCR duplicates are removed using UMI-tools (Smith **et al.** 2017).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

* python (developed on python 2.7.8)
* UMI-tools ()
* STAR (developed using STAR 2.5.0b)
* pandas
* HTSeq-count
* samtools

Additionally, an indexed reference genome in STAR format is required as well as the GTF that was used to create it.

### Installing

Currentlly, just run the script and ensure the dependencies above are satisfied.

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
	
