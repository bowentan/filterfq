# filterfq

The program for cleaning/filtereing raw fastq files, providing the following functionalities.

1. Automatically check quality system used in the raw fastq file(s)
2. Filter out reads that have a number of 'N' bases, low average quality or a number of low quality base
3. Convert quality system to specified system
4. Output statistical information of the raw and clean fastq reads, including distribution of read length, base, base quality
5. Multithread supported (up to 8)

## Getting Started

To use the program, you can fork a copy to your local machine by typing the following command in your work directory:

```
git clone https://github.com/bowentan/filterfq.git
```
or downloading the source codes or binary executables from the release page.

This program depends on Boost library of C++ with version C++11 or above, please make sure that your compiling environment satisfies those requirements for successful compilation.

For those whose system satisfies requirements and who download the source codes, use the following command to compile the program after unzipping it and entering the directory.

```
./configure --prefix=/your/package/path
make
make install
```
If your system cannot compile the source, please download the executable from the [Release](https://github.com/bowentan/filterfq/releases) page and tell us what problem you are facing in compiling so that we can fix it as soon as possible.

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Current versioning

v1.2.0

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/bowentan/filterfq/releases). 

## Authors

* **Bowen Tan** - *Initial work & maintainance* - [bowentan](https://github.com/bowentan)

## Acknowledgments

* Appreciate Wenlong Jia for the underlying mechanism of this program
* Appreciate Chang Xu for providing test data
