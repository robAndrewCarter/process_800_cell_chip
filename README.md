# Process the 800 cell IFC from Fluidigm

The fluidigm 800-cell IFC for the C1 has both barcodes and UMIs in read 1, whereas read 2 contains the actual sequence of a polyadenylated RNA. This script parses read1, generating cell-specific fastq files for read2 with UMIs appended to the read headers in read2. These Each demultiplexed sample in then aligned to a reference genome using STAR () and PCR duplicates are removed using UMI-tools ().

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

* python (developed on python 2.7.8)
* UMI-tools ()
* STAR (developed using STAR 2.5.0b)
* pandas
* HTSeq-count
* samtools

Additionally, an indexed reference genome in STAR format is required

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

The test folder contains a STAR-indexed genome of chromosomes - and - under the folder test/STAR_GENOME
Reads in the fastq folder test/FASTQ/ derive from 2 samples, each with coverage across two genes on each chromosome.
The test scripts run_tests.py under test/ runs the full script on this dataset and ensures that the resulting count matrix
assigns the reads to the correct samples, collapses reads using the UMIs, and derives the correct count matrix.


### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Rob Carter** - *Initial work* - (https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to - - for the mapping between cell barcode 6mer and row number

## References




