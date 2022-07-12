# coding: utf-8

# In[1]:

from Bio import SeqIO, SeqRecord, Seq
import pandas as pd
import regex, os, sys, re, subprocess, argparse, logging, shutil, pysam


# In[9]:


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""
    Generate sample-specific bam files from paired fastq files with cell-barcodes and UMIs."""
    )
    parser.add_argument("read1", type=str, help="path to read1")
    parser.add_argument("read2", type=str, help="path to read2")
    parser.add_argument("STAR_path", type=str, help="path to STAR index directory")
    parser.add_argument(
        "--sample_name",
        default="sample",
        type=str,
        help="name of sample. Used for file and directory naming",
    )
    parser.add_argument(
        "--barcode_mm",
        default=2,
        type=int,
        help="Number of mismatches allowed between cell barcodes",
    )
    parser.add_argument(
        "--barcode_start",
        default=0,
        type=int,
        help="start coordinate of barcode in read1",
    )
    parser.add_argument(
        "--barcode_end", default=5, type=int, help="end coordinate of barcode in read2"
    )
    parser.add_argument(
        "--umi_start", default=6, type=int, help="start coordinate of umi in read1"
    )
    parser.add_argument(
        "--umi_end", default=10, type=int, help="end coordinate of umi in read2"
    )
    parser.add_argument(
        "--output_dir", default="out", type=str, help="output directory"
    )
    parser.add_argument(
        "--output_format_string",
        default="foo{pos}moho",
        type=str,
        help="format string for naming files after splitting by barcode",
    )
    parser.add_argument(
        "--gff_filename",
        type=str,
        help="path to gff file used in STAR alignment and in htseq-count",
    )
    parser.add_argument(
        "--overwrite_option",
        default="exit",
        type=str,
        choices=["exit", "write", "read"],
        help="Behavior of program when previously created output files are encountered. Options = [exit (default), write, read]",
    )
    parser.add_argument(
        "--intervals_to_geneid_filename",
        type=str,
        default=None,
        help="BED formated file mapping exonic intervals to gene_ids. If None, the file is created from the gff_filename",
    )

    args = parser.parse_args()
    args.temp_dir = os.path.join(args.output_dir, args.sample_name)
    return args


def check_paths_and_setup(args):
    if os.path.exists(args.temp_dir):
        # logging.info("Temp directory {} exists. Removing it.".format(args.temp_dir))
        # try:
        #    shutil.rmtree(args.temp_dir)
        # except Exception:
        #    sys.exit("Cannot delete folder. Exiting")
        pass
    return True


# Set up global variables that would be given on command line if this is ever converted to a script

# In[130]:

barcode_to_row_dict = dict(
    CACGTA="ROW01",
    CTCACA="ROW02",
    TGCATC="ROW03",
    TCAGAC="ROW04",
    CGATGT="ROW05",
    TACTGC="ROW06",
    ATGCTC="ROW07",
    CATCTG="ROW08",
    GACTCA="ROW09",
    AGATCG="ROW10",
    ATCAGC="ROW11",
    GCTACA="ROW12",
    CAGATC="ROW13",
    CACAGT="ROW14",
    TACGAG="ROW15",
    CGACTA="ROW16",
    GCATCT="ROW17",
    AGCACT="ROW18",
    ACACTG="ROW19",
    CGTAGA="ROW20",
    ACTCGA="ROW21",
    ACATGC="ROW22",
    CGTCAT="ROW23",
    TGTACG="ROW24",
    GCAGTA="ROW25",
    TCACGT="ROW26",
    ACGTCA="ROW27",
    CTCGAT="ROW28",
    ATCGTG="ROW29",
    GCTGAT="ROW30",
    GTCTAC="ROW31",
    CATGCT="ROW32",
    TAGCAC="ROW33",
    GTGCAT="ROW34",
    TAGTCG="ROW35",
    TCTCAG="ROW36",
    CTAGTC="ROW37",
    TCTAGC="ROW38",
    ATGACG="ROW39",
    GAGCTA="ROW40",
    NONE="undetermined",
)


# In[17]:


def add_umis_to_fastq_reads(
    read1_fastq_filename, read2_fastq_filename, umi_start, umi_end, output_dir
):
    """
    Uses UMI-tools to add UMIs to read headers in both reads.
    The UMI sequences are assumed to be on read1.

    Returns:

    The paths to the forward and reverse reads of the umi-added fastq files.
    """
    r1_outfilename = os.path.join(
        output_dir,
        re.sub(
            ".fastq$", "_umi_labelled.fastq", os.path.basename(read1_fastq_filename)
        ),
    )
    r2_outfilename = os.path.join(
        output_dir,
        re.sub(
            ".fastq$", "_umi_labelled.fastq", os.path.basename(read2_fastq_filename)
        ),
    )
    if use_existing_files([r1_outfilename, r2_outfilename], args.overwrite_option):
        return (r1_outfilename, r2_outfilename)
    # Zero-based coordinates
    ssh_command = [
        "umi_tools",
        "extract",
        "--bc-pattern={}".format("X" * umi_start + "N" * (umi_end - umi_start + 1)),
        "-I",
        read1_fastq_filename,
        "-S",
        r1_outfilename,
        "--read2-in={}".format(read2_fastq_filename),
        "--read2-out={}".format(r2_outfilename),
    ]
    if subprocess.check_call(ssh_command):
        sys.exit("Problem running umi_tools")
    return (r1_outfilename, r2_outfilename)


def use_existing_files(file_list, write_option):
    if all([os.path.exists(_file) for _file in file_list]):
        if write_option == "exit":
            sys.exit(
                "At least one of {} exists. To overwrite, use overwrite_option = write. To use existing files, use overwrite_option = read".format(
                    ",".join(file_list)
                )
            )
        if write_option == "read":
            return True
        if write_option == "write":
            return False
    else:
        return False


# In[2]:


class FileManager:
    def __init__(self, mapping_dict, filename_format, overwrite_option="exit"):
        self.mapping_dict = mapping_dict
        self.filename_format = filename_format
        self.barcode_to_filehandle_dict = {}
        self.valid_barcodes = mapping_dict.keys()
        self.file_writing_stats_dict = {_barcode: 0 for _barcode in mapping_dict.keys()}
        self.num_files_written = 0
        possible_filenames = [
            self.filename_format.format(pos=_bc)
            for _bc in self.barcode_to_filehandle_dict.keys()
        ]
        self.nonempty_filenames = [
            _fn for _fn in possible_filenames if os.path.exists(_fn)
        ]
        self.status = None
        if len(self.nonempty_filenames) > 0:
            if overwrite_option == "exit":
                sys.exit(
                    "The following files already exist:{}. Use overwrite_option = read or overwrite_option = write".format(
                        ",".join(self.nonempty_filenames)
                    )
                )
            if overwrite_option == "read":
                self.status = "read"
            if overwrite_option == "write":
                self.nonempty_filenames = []
                self.status = "write"

    def write_seq(self, seq_rec, barcode):
        if not barcode in self.barcode_to_filehandle_dict:
            try:
                seq_filename = self.filename_format.format(
                    pos=self.mapping_dict[barcode]
                )
                self.nonempty_filenames.append(seq_filename)
                self.barcode_to_filehandle_dict[barcode] = open(seq_filename, mode="w")
            except Exception as e:
                sys.exit(
                    "Cannot create file {}".format(
                        self.filename_format.format(pos=self.mapping_dict[barcode])
                    )
                )
        try:
            SeqIO.write(
                format="fastq",
                handle=self.barcode_to_filehandle_dict[barcode],
                sequences=seq_rec,
            )
            self.file_writing_stats_dict[barcode] += 1
            self.num_files_written += 1
            if (self.num_files_written % 10000) == 0:
                logging.info(f"{self.num_files_written} files written")
        except Exception as e:
            temp_fh = self.barcode_to_filehandle_dict[barcode]
            sys.exit("Error writing to file {}".format(temp_fh.name))

    def get_stats(self):
        return self.file_writing_stats_dict

    def get_nonempty_filenames(self):
        return self.nonempty_filenames

    def close_filehandles(self):
        for _fh in self.barcode_to_filehandle_dict.values():
            _fh.close()


# In[3]:


def demultiplex_umi_labelled_fastq_files(
    cell_barcoded_umi_labelled_fastq_filename,
    umi_labelled_fastq_transcript_filename,
    mapping_dict,
    barcode_start,
    barcode_end,
    output_file_format_string,
    output_dir,
):
    file_manager = FileManager(
        filename_format=os.path.join(output_dir, output_file_format_string),
        mapping_dict=mapping_dict,
    )
    if file_manager.status == "read":
        return {
            "stats": None,
            "demultiplexed_umi_labelled_filenames": file_manager.get_nonempty_filenames(),
        }

    r1_gen = SeqIO.parse(cell_barcoded_umi_labelled_fastq_filename, "fastq")
    r2_gen = SeqIO.parse(umi_labelled_fastq_transcript_filename, "fastq")
    print(r2_gen)

    for (_barcoded_seqrec, _transcript_seqrec) in zip(r1_gen, r2_gen):
        barcode_seq = str(_barcoded_seqrec.seq[barcode_start : (barcode_end + 1)])
        matching_cell_barcode = []
        for _cell_barcode in barcode_to_row_dict.keys():
            m = regex.findall("(" + barcode_seq + "){s<=1}", _cell_barcode)
            if m:
                matching_cell_barcode.append(_cell_barcode)
        num_matches = len(matching_cell_barcode)
        if num_matches == 1:
            _transcript_seqrec.id = ":".join(
                [
                    "BCM",
                    matching_cell_barcode[0],
                    "BC",
                    barcode_seq,
                    _transcript_seqrec.id,
                ]
            )
            file_manager.write_seq(_transcript_seqrec, matching_cell_barcode[0])
        elif num_matches == 0:
            _transcript_seqrec.id = ":".join(
                ["BCM", "NONE", "BC", barcode_seq, _transcript_seqrec.id]
            )
            file_manager.write_seq(_transcript_seqrec, "NONE")
        elif num_matches > 1:
            sys.exit("sequence matches more than two barcodes")
    file_manager.close_filehandles()
    return {
        "stats": file_manager.get_stats(),
        "demultiplexed_umi_labelled_filenames": file_manager.get_nonempty_filenames(),
    }


# In[ ]:


def align_sequences(fastq_filenames_list, star_path, output_dir):
    """
    Aligns a list of fastq files to the specified genomes using STAR's 2pass alignment.
    The resulting BAM files are also indexed

    Returns:

    a list of bam file paths corresponding to the 2nd pass of the STAR alignment
    """

    aligned_star_pass2_bam_filenames_list = []
    for _filename in fastq_filenames_list:
        _filename_prefix = re.sub("\..+", "", os.path.basename(_filename))
        _temp_star1_dir = os.path.join(output_dir, _filename_prefix + "_star1")
        _temp_star2_dir = os.path.join(output_dir, _filename_prefix + "_star2")
        _filename_star1_prefix = os.path.join(_temp_star1_dir, _filename_prefix)
        _filename_star2_prefix = os.path.join(_temp_star2_dir, _filename_prefix)
        _star_pass2_bamfile = _filename_star2_prefix + "Aligned.sortedByCoord.out.bam"
        if not use_existing_files([_temp_star1_dir], args.overwrite_option):
            subprocess.check_call(["mkdir", "-p", _temp_star1_dir])
            # subprocess.check_call(["cd", _temp_dir])
            # subprocess.check_call(['module','load','star/2.5.0b'])
            subprocess.check_call(
                [
                    "STAR",
                    "--runThreadN",
                    "1",
                    "--runMode",
                    "alignReads",
                    "--genomeDir",
                    star_path,
                    "--readFilesIn",
                    _filename,
                    "--outFileNamePrefix",
                    _filename_star1_prefix,
                    "--outSAMtype",
                    "BAM",
                    "SortedByCoordinate",
                    "--quantMode",
                    "GeneCounts",
                    "--outSAMstrandField",
                    "intronMotif",
                ]
            )
        if not use_existing_files([_temp_star2_dir], args.overwrite_option):
            out_tab_filename = _filename_star1_prefix + "SJ.out.tab"
            subprocess.check_call(["mkdir", "-p", _temp_star2_dir])
            if os.path.exists(out_tab_filename):
                subprocess.check_call(
                    [
                        "STAR",
                        "--runThreadN",
                        "1",
                        "--runMode",
                        "alignReads",
                        "--genomeDir",
                        star_path,
                        "--readFilesIn",
                        _filename,
                        "--outFileNamePrefix",
                        _filename_star2_prefix,
                        "--outSAMtype",
                        "BAM",
                        "SortedByCoordinate",
                        "--quantMode",
                        "GeneCounts",
                        "--sjdbFileChrStartEnd",
                        out_tab_filename,
                        "--outSAMstrandField",
                        "intronMotif",
                    ]
                )
            else:
                sys.exit("Could not find tab file for STAR pass 2")
        if os.path.exists(_star_pass2_bamfile):
            aligned_star_pass2_bam_filenames_list.append(_star_pass2_bamfile)
            try:
                subprocess.check_call(["samtools", "index", _star_pass2_bamfile])
            except Exception:
                sys.exit("Could not run samtools index {}".format(_star_pass2_bamfile))
        else:
            sys.exit("Could not find bam file {}".format(_star_pass2_bamfile))
    return aligned_star_pass2_bam_filenames_list


def tag_bams(aligned_bam_filenames, intervals_to_geneid_filename, output_dir):
    tagged_bam_filenames_list = []
    for _bam_filename in aligned_bam_filenames:
        tagged_bam_filename = os.path.join(
            output_dir,
            re.sub("\.[^.]+$", "_tagged.bam", os.path.basename(_bam_filename)),
        )
        run_command = [
            "bedtools",
            "tag",
            "-names",
            "-tag",
            "GN",
            "-i",
            _bam_filename,
            "-files",
            intervals_to_geneid_filename,
        ]
        try:
            with open(tagged_bam_filename, "w") as outfile:
                subprocess.call(run_command, stdout=outfile)
        except:
            sys.exit("Cannot run {}".format(run_command))
        temp_tagged_bam_filename = tagged_bam_filename + ".tmp"
        try:
            in_bam_obj = pysam.AlignmentFile(tagged_bam_filename, mode="rb")
            out_bam_obj = pysam.AlignmentFile(
                temp_tagged_bam_filename, template=in_bam_obj, mode="wb"
            )
        except:
            sys.exit("Can't process bams")
        for _read in in_bam_obj:
            if _read.has_tag("GN"):
                _read.set_tag("GN", ",".join(list(set(_read.get_tag("GN").split(",")))))
                out_bam_obj.write(_read)
        in_bam_obj.close()
        out_bam_obj.close()
        try:
            run_command = ["mv", temp_tagged_bam_filename, tagged_bam_filename]
            subprocess.check_call(run_command)
        except:
            sys.exit("Cannot mv file")
        sort_and_index_bam_filename(tagged_bam_filename)
        tagged_bam_filenames_list.append(tagged_bam_filename)
    return tagged_bam_filenames_list


def sort_and_index_bam_filename(bam_filename):
    temp_sorted_filename = (
        re.sub("\.[^.]+$", "", os.path.basename(bam_filename)) + "_sorted.bam.tmp"
    )
    try:
        subprocess.check_call(
            ["samtools", "sort", "-o", temp_sorted_filename, bam_filename]
        )
    except Exception:
        sys.exit("Failed to sort {}".format(bam_filename))
    try:
        subprocess.check_call(["mv", temp_sorted_filename, bam_filename])
    except:
        sys.exit("Failed to rename {}".format(temp_sorted_filename))
    try:
        subprocess.check_call(["samtools", "index", bam_filename])
    except:
        sys.exit("Failed to index {}".format(bam_filename))
    return True


def dedup_bam_files(bam_tagged_filenames_list, output_dir):
    deduped_bam_filenames_list = []
    for _bam_filename in bam_tagged_filenames_list:
        _deduped_bam_filename = os.path.join(
            output_dir,
            re.sub("\.bam$", "_deduped.bam", os.path.basename(_bam_filename)),
        )
        _deduped_log_filename = os.path.join(
            output_dir,
            re.sub("\.bam$", "_deduped.log", os.path.basename(_bam_filename)),
        )
        try:
            subprocess.check_call(
                [
                    "umi_tools",
                    "dedup",
                    "-v",
                    "10",
                    "--gene-tag=GN",
                    "--per-gene",
                    "-I",
                    _bam_filename,
                    "-S",
                    _deduped_bam_filename,
                    "-L",
                    _deduped_log_filename,
                ]
            )
            deduped_bam_filenames_list.append(_deduped_bam_filename)
        except Exception:
            sys.exit("Unable to successfully run umi_tools dedup")
    return deduped_bam_filenames_list


def get_gene_counts_from_tagged_bam_files(
    bam_filenames_list, ordered_genes_list, output_dir
):
    counts_df_list = []
    for _bam_filename in bam_filenames_list:
        try:
            bam_obj = pysam.AlignmentFile(_bam_filename, "rb")
            read_to_gene_dict = {}
            for _read in bam_obj:
                if _read.query_name in read_to_gene_dict:
                    read_to_gene_dict[_read.query_name].append(_read.get_tag("GN"))
                else:
                    read_to_gene_dict[_read.query_name] = [_read.get_tag("GN")]
        except:
            sys.exit("problem reading bam file {}".format(_bam_filename))

        gene_count_dict = {_gene: 0 for _gene in ordered_genes_list}

        singlet_count = 0
        multiplet_count = 0

        for _read_id, gene_list in read_to_gene_dict.iteritems():
            if len(gene_list) > 1:
                multiplet_count += 1
            elif len(gene_list) == 1:
                singlet_count += 1
                gene_count_dict[gene_list[0]] += 1

        counts_df = (
            pd.DataFrame.from_dict(gene_count_dict, orient="index")
            .sort_index()
            .rename(
                columns={
                    0: os.path.join(
                        output_dir, re.sub("\..+", "", os.path.basename(_bam_filename))
                    )
                }
            )
        )
        counts_df_list.append(counts_df)
        try:
            counts_df.to_csv(
                os.path.join(
                    output_dir,
                    re.sub("\..+", "", os.path.basename(_bam_filename))
                    + "_gene_counts.tsv",
                ),
                sep="\t",
                header=False,
                index=False,
            )
        except:
            sys.exit("Cannot open file for writing")
    return counts_df_list


def merge_and_write_counts(count_df_list, sample_name, output_dir):
    # merge count DataFrames and write to file as a matrix
    merged_count_df = pd.concat(counts_df_list, axis=1, ignore_index=False).T
    try:
        merged_count_df.to_csv(
            os.path.join(output_dir, sample_name + "_samples_gene_counts.tsv"), sep="\t"
        )
    except Exception:
        sys.exit("Cannot write matrix file")
    return True


def get_intervals_to_geneid_filename(
    intervals_to_geneid_filename, gff_filename, output_dir
):
    if not intervals_to_geneid_filename:
        try:
            temp_exons_df = pd.read_csv(gff_filename, sep="\t", header=None)
        except:
            sys.exit("Cannot open {}".format(gff_filename))
        temp_exons_df = temp_exons_df[
            temp_exons_df[8].str.contains("gene_id") & (temp_exons_df[2] == "exon")
        ]
        temp_exons_df.loc[:, "gene_id"] = temp_exons_df[8].str.extract(
            'gene_id "(?P<gene_id>[^"]+)"', expand=True
        )
        intervals_bed_filename = os.path.join(
            output_dir,
            re.sub("\..+$", "_exons_to_geneid.bed", os.path.basename(gff_filename)),
        )
        try:
            temp_exons_df.loc[:, [0, 3, 4, "gene_id"]].to_csv(
                intervals_bed_filename, sep="\t", header=False, index=False
            )
        except:
            sys.exit("Cannot write {}".format())
        return intervals_bed_filename
    else:
        return intervals_to_geneid_filename


def get_gene_counts_using_htseq(bam_filenames_list, gff_filename, output_dir):
    htseq_filenames_list = []
    for _bam_filename in bam_filenames_list:
        _htseq_filename = os.path.join(
            output_dir, re.sub("\.bam$", "_htseq.tsv", os.path.basename(_bam_filename))
        )
        # _htseq_log_filename = os.path.join(output_dir, re.sub("\.bam$", "_htseq.log", os.path.basename(_bam_filename)))
        _command = ["htseq-count", "-f", "bam", "-s", "no", _bam_filename, gff_filename]
        try:
            _htseq_output = subprocess.check_output(_command)
            print(_htseq_output, _htseq_filename)
        except Exception:
            sys.exit(
                "Unable to successfully run htseq-count: {}".format(" ".join(_command))
            )
        try:
            _fh = open(_htseq_filename, "w")
            _fh.write(_htseq_output.decode("utf-8"))
            _fh.close()
        except:
            sys.exit("Unable to successfully write htseq files")
        htseq_filenames_list.append(_htseq_filename)
    return htseq_filenames_list


def merge_and_write_htseq_counts(htseq_filenames_list, output_dir):
    df_list = []
    for _htseq_filename in htseq_filenames_list:
        htseq_counts_df = pd.read_csv(
            _htseq_filename,
            sep="\t",
            index_col=0,
            header=None,
            names=[re.sub("\..+$", "", os.path.basename(_htseq_filename))],
        )
        htseq_counts_df.index = [
            re.sub("\..+$", "", _entry) for _entry in htseq_counts_df.index
        ]
        df_list.append(htseq_counts_df)
    merged_df = pd.concat(df_list, axis=1)
    merged_df.T.to_csv(
        os.path.join(output_dir, args.sample_name + "_samples_gene_counts.tsv"),
        sep="\t",
    )


# In[ ]:

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    args = parse_arguments()
    check_paths_and_setup(args)
    # create temp directory
    try:
        if not os.path.exists(args.temp_dir):
            os.makedirs(args.temp_dir)
    except Exception:
        sys.exit(f"Unable to create directory {args.temp_dir}")
    (r1_umi_fastq_path, r2_umi_fastq_path) = add_umis_to_fastq_reads(
        args.read1, args.read2, args.umi_start, args.umi_end, args.temp_dir
    )
    logging.info(
        "UMI sequences added to files {} and {} from {} and {}, respectively".format(
            r1_umi_fastq_path, r2_umi_fastq_path, args.read1, args.read2
        )
    )

    file_manager_info_dict = demultiplex_umi_labelled_fastq_files(
        r1_umi_fastq_path,
        r2_umi_fastq_path,
        barcode_to_row_dict,
        args.barcode_start,
        args.barcode_end,
        args.output_format_string,
        args.temp_dir,
    )
    aligned_bam_filenames = align_sequences(
        file_manager_info_dict["demultiplexed_umi_labelled_filenames"],
        args.STAR_path,
        args.temp_dir,
    )
    intervals_to_geneid_filename = get_intervals_to_geneid_filename(
        args.intervals_to_geneid_filename, args.gff_filename, args.temp_dir
    )
    tagged_bam_filenames = tag_bams(
        aligned_bam_filenames, intervals_to_geneid_filename, args.temp_dir
    )
    deduped_bam_filenames = dedup_bam_files(tagged_bam_filenames, args.temp_dir)
    htseq_count_filenames = get_gene_counts_using_htseq(
        deduped_bam_filenames, args.gff_filename, args.temp_dir
    )
    merge_and_write_htseq_counts(htseq_count_filenames, args.output_dir)
