# ClipKIT: The **Clip**ping and **K**eeping **I**nformation **T**rimmer

ClipKIT is a fast and flexible alignment trimming tool. ClipKIT keeps phylogenetically informative sites while removing others
that display characteristics of poor phylogenetic signal.<br />
Add journal and citation information
<br /><br />


---


## Guide
[Installation](#installation)<br />
[Quick Start](#quick-start)<br />
[Detailed Usage](#detailed-usage)


---


### Installation
To install, use the following commands:
```
```
<br />
To install from source, use the following commands:
```
```
<br />

### Quick Start
ClipKIT requires only two arguments, the input and output. To use ClipKIT in its simpliest form, use the following command:
```
clipkit -i $input -o $output
```
where $input is your input alignment file and $output is the name of your output file

<br />

### Detailed Usage
This section described the various features and options of ClipKIT.
[Modes](#modes)<br />
[Log](#log)<br />
[Complementary](#complementary)<br />
[Miscellaneous options](#miscellaneous-options)

<br />

#### Modes
ClipKIT can run with three different modes (kpi, gappy, kpi-gappy), which are specified with the -m/--mode argument.<br /> 
*Default: 'gappy'*<br />
* kpi will trim all sites that are not parsimony informative<br />
* gappy will remove all sites that are above a threshold of gappyness (default: 0.9)<br />
* kpi-gappy is the combination of kpi- and gappy-based trimming<br />
```
# kpi-based trimming
clipkit -i $input -o $output -m kpi

# gappy-based trimming
clipkit -i $input -o $output -m gappy

# kpi- and gappy-based trimming
clipkit -i $input -o $output -m kpi-gappy 
```

<br />

#### Log
It can be very useful to have information about the each position in an alignment. For example, this information could be used in alignment diagnostics, fine-tuning of trimming parameters, etc. To create the log file, use the -l/--log option. Using the -l/--log option will create a four column file with the suffix '.log'. 
* col1: position in the alignment (starting at 1)
* col2: reports if site was trimmed or kept (trim or keep, respectively)
* col3: reports if the site is parsimony informative or not (PI or nPI, respectively)
* col4: reports the gappyness of the position (number of gaps / entries in alignment)
*Default: off*<br />
```
clipkit -i $input -o $output -l
```
This will result in an additional output file named $output.log

<br />

### Complementary
Having an alignment of the sequences that were trimmed can be useful for other analyses. To obtain an alignment of the sequences that were trimmed, use the -c/--complementary option. Using the -c/--complementary option will create a file with the suffix '.complement'.
*Default: off*<br />

```
clipkit -i $input -o $output -c
```
This will result in an additional output file named $output.complement

<br />

#### Miscellaneous options
| Option        | Usage and meaning |
| ------------- | ------------------ |
| -h/--help     | Print help message |
| -v/--version  | Print software version |
| -g/--gaps     | Specify gappyness threshold (between 0 and 1). *Default: 0.9* |
| -if/--input_file_format | Specify input file format. Accepted file formats are: fasta, clustal, maf, mauve, phylip, phylip-sequential, phylip-relaxed, stockholm. *Default: auto-detect* |
| -of/--input_file_format | Specify output file format. Accepted file formats are: fasta, clustal, maf, mauve, phylip, phylip-sequential, phylip-relaxed, stockholm. *Default: input file type* |

<br />

##### Accepted file formats formats
[fasta](https://en.wikipedia.org/wiki/FASTA_format), [clustal](http://meme-suite.org/doc/clustalw-format.html), [maf](http://www.bx.psu.edu/~dcking/man/maf.xhtml), [mauve](http://darlinglab.org/mauve/user-guide/files.html), [phylip](http://scikit-bio.org/docs/0.2.3/generated/skbio.io.phylip.html), [phylip-sequential](http://rosalind.info/glossary/phylip-format/), [phylip-relaxed](https://www.hiv.lanl.gov/content/sequence/FORMAT_CONVERSION/FormatExplain.html), [stockholm](https://en.wikipedia.org/wiki/Stockholm_format)

<br />

## Authors
* [Jacob Steenwyk](https://jlsteenwyk.github.io/)<br />
* [Thomas Buida](www.tjbiii.com)<br />
* Others
<br />

## Developers
* [Jacob Steenwyk](https://jlsteenwyk.github.io/)<br />
* [Thomas Buida](www.tjbiii.com)<br />
<br />


